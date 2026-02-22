"""
Main window for Dorothy application.
"""

import os
from pathlib import Path

from PyQt6.QtWidgets import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLineEdit,
    QLabel,
    QStackedWidget,
    QScrollArea,
    QFrame,
    QComboBox,
    QSpinBox,
    QCheckBox,
    QGroupBox,
    QProgressBar,
    QFileDialog,
    QMessageBox,
    QDialog,
    QDialogButtonBox,
)
from PyQt6.QtCore import Qt, QCoreApplication, QThread, pyqtSignal
from PyQt6.QtGui import QFont, QDesktopServices, QPixmap, QIcon
from PyQt6.QtCore import QUrl

from dorothy.core.cod_search import CODSearch, MoleculeResult
from dorothy.core.cif_parser import MoleculeStructure
from dorothy.core.generator import GenerationPipeline, GenerationSettings, GenerationResult
from dorothy.core.xtb_manager import is_xtb_installed, download_xtb, get_download_url, get_xtb_install_instructions
from dorothy.core.density import create_density_cube_from_structure
from dorothy.ui.molecule_viewer import MoleculeViewer
from dorothy.ui.slice_explorer import SliceExplorer
from dorothy.ui import styles as S


# =============================================================================
# Font helper — uses the design system font at a given size/weight
# =============================================================================

def _font(size: int = S.SIZE_BODY, bold: bool = False) -> QFont:
    f = QFont(S.FONT_FAMILY, size)
    f.setBold(bold)
    return f


# =============================================================================
# Worker threads (unchanged logic)
# =============================================================================

class SearchWorker(QThread):
    """Background thread for COD search."""
    finished = pyqtSignal(list, bool)  # results, is_online

    def __init__(self, searcher: CODSearch, query: str):
        super().__init__()
        self.query = query
        self.searcher = searcher

    def run(self):
        results, is_online = self.searcher.search(self.query)
        self.finished.emit(results, is_online)


class GenerationWorker(QThread):
    """Background thread for PDF generation."""
    progress = pyqtSignal(str, int)  # message, percent
    finished = pyqtSignal(object)  # GenerationResult

    def __init__(self, structure: MoleculeStructure, output_dir: Path, settings: GenerationSettings):
        super().__init__()
        self.structure = structure
        self.output_dir = output_dir
        self.settings = settings

    def run(self):
        pipeline = GenerationPipeline(
            progress_callback=lambda msg, pct: self.progress.emit(msg, pct)
        )
        result = pipeline.generate(self.structure, self.output_dir, self.settings)
        self.finished.emit(result)


class XtbDownloadWorker(QThread):
    """Background thread for xTB download."""
    progress = pyqtSignal(int, int)  # downloaded, total
    finished = pyqtSignal(bool)  # success

    def run(self):
        success = download_xtb(
            progress_callback=lambda downloaded, total: self.progress.emit(downloaded, total)
        )
        self.finished.emit(success)


class XtbDensityWorker(QThread):
    """Background thread for xTB density calculation."""
    progress = pyqtSignal(str)  # status message
    finished = pyqtSignal(object, object, object)  # promolecule, molecular, deformation

    def __init__(self, structure: MoleculeStructure, plane_definition=None):
        super().__init__()
        self.structure = structure
        self.plane_definition = plane_definition

    def run(self):
        from dorothy.core.xtb_manager import is_xtb_installed, run_xtb_density
        from dorothy.core.density import (
            parse_cube_file,
            calculate_deformation_density,
            create_density_cube_from_structure
        )
        import tempfile

        promolecule = None
        molecular = None
        deformation = None

        try:
            self.progress.emit("Generating promolecule density...")
            promolecule = create_density_cube_from_structure(
                self.structure,
                resolution="coarse",
                align_to_principal_axes=(self.plane_definition is None),
                plane_definition=self.plane_definition
            )

            if is_xtb_installed():
                self.progress.emit("Running xTB calculation...")
                with tempfile.TemporaryDirectory() as tmp_dir:
                    cube_path = run_xtb_density(
                        self.structure,
                        Path(tmp_dir),
                        lambda msg: self.progress.emit(msg),
                        plane_definition=self.plane_definition,
                    )

                    if cube_path:
                        self.progress.emit("Processing density data...")
                        molecular = parse_cube_file(cube_path)
                        if molecular:
                            promolecule, deformation = calculate_deformation_density(
                                molecular,
                                self.structure
                            )

            self.progress.emit("Ready")

        except Exception as e:
            self.progress.emit(f"Error: {e}")

        self.finished.emit(promolecule, molecular, deformation)


# =============================================================================
# xTB Download Dialog
# =============================================================================

class XtbDownloadDialog(QDialog):
    """Dialog for downloading/installing xTB."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Install xTB")
        self.setMinimumWidth(450)
        self.setModal(True)

        self.download_worker = None
        self.can_auto_download = get_download_url() is not None

        layout = QVBoxLayout(self)
        layout.setSpacing(15)

        title = QLabel("xTB Not Installed")
        title.setFont(_font(S.SIZE_SUBTITLE, bold=True))
        title.setStyleSheet(S.label_title())
        layout.addWidget(title)

        explanation = QLabel(
            "xTB is required for deformation density calculation.\n\n"
            "Without xTB, only promolecule density (atomic positions) "
            "will be generated. With xTB, you also get deformation density "
            "showing chemical bonding."
        )
        explanation.setWordWrap(True)
        explanation.setFont(_font())
        explanation.setStyleSheet(S.label_secondary())
        layout.addWidget(explanation)

        instructions = get_xtb_install_instructions()
        install_label = QLabel(instructions)
        install_label.setWordWrap(True)
        install_label.setStyleSheet(
            f"background-color: {S.BG_CARD}; padding: 10px; "
            f"font-family: monospace; border-radius: 6px; color: {S.TEXT};"
        )
        install_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(install_label)

        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setStyleSheet(S.progress_bar())
        self.progress_bar.hide()
        layout.addWidget(self.progress_bar)

        self.status_label = QLabel()
        self.status_label.setFont(_font())
        self.status_label.setStyleSheet(S.label_secondary())
        self.status_label.hide()
        layout.addWidget(self.status_label)

        button_layout = QHBoxLayout()

        self.skip_btn = QPushButton("Skip (Promolecule Only)")
        self.skip_btn.setFont(_font())
        self.skip_btn.setStyleSheet(S.btn_secondary())
        self.skip_btn.clicked.connect(self.reject)
        button_layout.addWidget(self.skip_btn)

        if self.can_auto_download:
            self.download_btn = QPushButton("Download xTB")
            self.download_btn.setFont(_font(bold=True))
            self.download_btn.setStyleSheet(S.btn_primary())
            self.download_btn.clicked.connect(self._start_download)
            button_layout.addWidget(self.download_btn)
        else:
            self.check_btn = QPushButton("I've Installed It")
            self.check_btn.setFont(_font(bold=True))
            self.check_btn.setStyleSheet(S.btn_primary())
            self.check_btn.clicked.connect(self._check_installation)
            button_layout.addWidget(self.check_btn)

        layout.addLayout(button_layout)

    def _check_installation(self):
        if is_xtb_installed():
            self.accept()
        else:
            self.status_label.show()
            self.status_label.setText("xTB not found. Please install it and try again.")
            self.status_label.setStyleSheet(S.label_error())

    def _start_download(self):
        self.download_btn.setEnabled(False)
        self.skip_btn.setEnabled(False)
        self.progress_bar.show()
        self.progress_bar.setValue(0)
        self.status_label.show()
        self.status_label.setText("Downloading...")

        self.download_worker = XtbDownloadWorker()
        self.download_worker.progress.connect(self._on_progress)
        self.download_worker.finished.connect(self._on_finished)
        self.download_worker.start()

    def _on_progress(self, downloaded: int, total: int):
        if total > 0:
            percent = int(100 * downloaded / total)
            self.progress_bar.setValue(percent)
            mb_downloaded = downloaded / (1024 * 1024)
            mb_total = total / (1024 * 1024)
            self.status_label.setText(f"Downloading... {mb_downloaded:.1f} / {mb_total:.1f} MB")

    def _on_finished(self, success: bool):
        if success:
            self.status_label.setText("Download complete!")
            self.status_label.setStyleSheet(S.label_success())
            self.accept()
        else:
            self.status_label.setText("Download failed. Please install manually.")
            self.status_label.setStyleSheet(S.label_error())
            self.download_btn.setEnabled(True)
            self.skip_btn.setEnabled(True)


# =============================================================================
# Main Window
# =============================================================================

class MainWindow(QMainWindow):
    """Main application window for Dorothy."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Dorothy")
        self.setMinimumSize(800, 600)
        icon_path = Path(__file__).parent.parent.parent / "logo" / "dorothy_logo_icon.png"
        if icon_path.exists():
            self.setWindowIcon(QIcon(str(icon_path)))
        self.showMaximized()

        self.searcher = CODSearch()
        self.search_worker: SearchWorker | None = None
        self.generation_worker: GenerationWorker | None = None
        self.xtb_density_worker: XtbDensityWorker | None = None
        self.current_results = []
        self.selected_molecule: MoleculeResult | None = None
        self.selected_structure: MoleculeStructure | None = None
        self.last_result: GenerationResult | None = None
        self.selection_manager = None

        central_widget = QWidget()
        central_widget.setObjectName("dorothy_central")
        central_widget.setStyleSheet("#dorothy_central { background-color: #ffffff; }")
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        self.stacked_widget = QStackedWidget()
        main_layout.addWidget(self.stacked_widget)

        self.home_screen = self._create_home_screen()
        self.results_screen = self._create_results_screen()
        self.preview_screen = self._create_preview_screen()
        self.processing_screen = self._create_processing_screen()
        self.complete_screen = self._create_complete_screen()
        self.stacked_widget.addWidget(self.home_screen)
        self.stacked_widget.addWidget(self.results_screen)
        self.stacked_widget.addWidget(self.preview_screen)
        self.stacked_widget.addWidget(self.processing_screen)
        self.stacked_widget.addWidget(self.complete_screen)

        self.stacked_widget.setCurrentWidget(self.home_screen)

    # =========================================================================
    # Home Screen
    # =========================================================================

    def _create_home_screen(self) -> QWidget:
        screen = QWidget()
        screen.setObjectName("home_screen")
        screen.setStyleSheet("#home_screen { background-color: #ffffff; }")
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setSpacing(0)

        layout.addStretch(2)

        # Logo
        logo_label = QLabel()
        logo_path = Path(__file__).parent.parent.parent / "logo" / "dorothy_logo.png"
        if logo_path.exists():
            pixmap = QPixmap(str(logo_path))
            scaled = pixmap.scaledToWidth(350, Qt.TransformationMode.SmoothTransformation)
            logo_label.setPixmap(scaled)
        else:
            logo_label.setText(self.tr("Dorothy"))
            logo_label.setFont(_font(32, bold=True))
        logo_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(logo_label)

        layout.addSpacing(8)

        # Subtitle
        subtitle = QLabel(self.tr("A tool for discovery and learning"))
        subtitle.setFont(_font(16))
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        subtitle.setStyleSheet(f"color: {S.TEXT_SECONDARY};")
        layout.addWidget(subtitle)

        layout.addSpacing(30)

        # Google-style search bar — single wide pill, centered
        search_row = QHBoxLayout()
        search_row.setContentsMargins(0, 0, 0, 0)
        search_row.addStretch(1)

        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText(self.tr("Search molecule (e.g. aspirin, benzene, caffeine)..."))
        self.search_input.setFont(_font(16))
        self.search_input.setStyleSheet(f"""
            QLineEdit {{
                background-color: white;
                color: {S.TEXT};
                border: 1px solid {S.BORDER};
                border-radius: 24px;
                padding: 12px 24px;
                font-size: 16pt;
                font-family: "{S.FONT_FAMILY}";
            }}
            QLineEdit:focus {{
                border-color: {S.ACCENT};
                border-width: 2px;
            }}
        """)
        self.search_input.setFixedWidth(580)
        self.search_input.returnPressed.connect(self._on_search)
        search_row.addWidget(self.search_input)

        self.search_button = QPushButton(self.tr("Search"))
        self.search_button.setObjectName("search_btn")
        self.search_button.setStyleSheet(f"""
            #search_btn {{
                background-color: {S.ACCENT};
                color: white;
                border: 2px solid {S.ACCENT};
                border-radius: 24px;
                padding: 12px 28px;
                font-size: 16pt;
                font-family: "{S.FONT_FAMILY}";
                font-weight: bold;
            }}
            #search_btn:hover {{
                background-color: {S.ACCENT_HOVER};
                border-color: {S.ACCENT_HOVER};
            }}
            #search_btn:pressed {{
                background-color: {S.ACCENT_PRESSED};
                border-color: {S.ACCENT_PRESSED};
            }}
        """)
        self.search_button.clicked.connect(self._on_search)
        search_row.addSpacing(10)
        search_row.addWidget(self.search_button)

        search_row.addStretch(1)
        layout.addLayout(search_row)

        layout.addSpacing(12)

        # Status label (shows "Searching..." etc.)
        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.status_label.setFont(_font(S.SIZE_SUBTITLE))
        self.status_label.setStyleSheet(f"color: {S.TEXT_SECONDARY};")
        layout.addWidget(self.status_label)

        layout.addStretch(3)

        # Footer
        footer = QLabel(
            self.tr("Honoring Dorothy Hodgkin's pioneering work in X-ray crystallography")
        )
        footer.setAlignment(Qt.AlignmentFlag.AlignCenter)
        footer.setFont(_font(13))
        footer.setStyleSheet(f"color: {S.TEXT_MUTED};")
        layout.addWidget(footer)

        return screen

    # =========================================================================
    # Results Screen
    # =========================================================================

    def _create_results_screen(self) -> QWidget:
        screen = QWidget()
        screen.setObjectName("results_screen")
        screen.setStyleSheet("#results_screen { background-color: #ffffff; }")
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        header = QHBoxLayout()
        back_btn = QPushButton(self.tr("< Back"))
        back_btn.setFont(_font())
        back_btn.setStyleSheet(S.btn_secondary())
        back_btn.clicked.connect(self._go_home)
        header.addWidget(back_btn)

        self.results_title = QLabel(self.tr("Results"))
        self.results_title.setFont(_font(S.SIZE_TITLE, bold=True))
        self.results_title.setStyleSheet(S.label_title())
        header.addWidget(self.results_title)
        header.addStretch()

        self.online_status = QLabel()
        self.online_status.setFont(_font(S.SIZE_SMALL))
        header.addWidget(self.online_status)

        layout.addLayout(header)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        self.results_container = QWidget()
        self.results_layout = QVBoxLayout(self.results_container)
        self.results_layout.setSpacing(10)
        self.results_layout.addStretch()

        scroll.setWidget(self.results_container)
        layout.addWidget(scroll)

        return screen

    # =========================================================================
    # Preview Screen
    # =========================================================================

    def _create_preview_screen(self) -> QWidget:
        screen = QWidget()
        screen.setObjectName("preview_screen")
        screen.setStyleSheet("#preview_screen { background-color: #ffffff; }")
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(15)

        # Header
        header = QHBoxLayout()
        back_btn = QPushButton(self.tr("< Back"))
        back_btn.setFont(_font())
        back_btn.setStyleSheet(S.btn_secondary())
        back_btn.clicked.connect(self._go_to_results)
        header.addWidget(back_btn)

        self.preview_title = QLabel(self.tr("Molecule"))
        self.preview_title.setFont(_font(S.SIZE_TITLE, bold=True))
        self.preview_title.setStyleSheet(S.label_title())
        header.addWidget(self.preview_title)
        header.addStretch()

        layout.addLayout(header)

        # Main content
        content = QHBoxLayout()
        viewer_container = QVBoxLayout()

        # View toggle (2D / 3D)
        view_toggle = QHBoxLayout()
        self.view_2d_btn = QPushButton(self.tr("2D Structure"))
        self.view_2d_btn.setFont(_font())
        self.view_2d_btn.setCheckable(True)
        self.view_2d_btn.setChecked(True)
        self.view_2d_btn.clicked.connect(lambda: self._switch_view("2d"))
        self.view_2d_btn.setStyleSheet(S.btn_toggle())
        view_toggle.addWidget(self.view_2d_btn)

        self.view_3d_btn = QPushButton(self.tr("3D Slices"))
        self.view_3d_btn.setFont(_font())
        self.view_3d_btn.setCheckable(True)
        self.view_3d_btn.clicked.connect(lambda: self._switch_view("3d"))
        self.view_3d_btn.setStyleSheet(S.btn_toggle())
        view_toggle.addWidget(self.view_3d_btn)

        view_toggle.addSpacing(20)

        self.clear_selection_btn = QPushButton(self.tr("Clear Selection"))
        self.clear_selection_btn.setFont(_font())
        self.clear_selection_btn.setStyleSheet(S.btn_secondary())
        self.clear_selection_btn.clicked.connect(self._clear_atom_selection)
        self.clear_selection_btn.setToolTip("Clear atom selection for plane definition")
        view_toggle.addWidget(self.clear_selection_btn)

        self.reset_plane_btn = QPushButton(self.tr("Reset Plane"))
        self.reset_plane_btn.setFont(_font())
        self.reset_plane_btn.setStyleSheet(S.btn_secondary())
        self.reset_plane_btn.clicked.connect(self._reset_slice_plane)
        self.reset_plane_btn.setToolTip("Reset to default slicing orientation")
        view_toggle.addWidget(self.reset_plane_btn)

        view_toggle.addStretch()
        viewer_container.addLayout(view_toggle)

        # Viewer stack (2D / 3D)
        self.viewer_stack = QStackedWidget()

        self.molecule_viewer = MoleculeViewer()
        self.molecule_viewer.setMinimumSize(350, 350)
        self.viewer_stack.addWidget(self.molecule_viewer)

        self.slice_explorer = SliceExplorer()
        self.slice_explorer.setMinimumSize(350, 350)
        self.viewer_stack.addWidget(self.slice_explorer)

        viewer_container.addWidget(self.viewer_stack)

        # Molecule info
        self.mol_info_label = QLabel()
        self.mol_info_label.setFont(_font())
        self.mol_info_label.setStyleSheet(S.label_secondary())
        self.mol_info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        viewer_container.addWidget(self.mol_info_label)

        content.addLayout(viewer_container, stretch=1)
        layout.addLayout(content)

        # Collapsible PDF Export section
        self.pdf_toggle_btn = QPushButton(self.tr("PDF Export Settings  >"))
        self.pdf_toggle_btn.setFont(_font(bold=True))
        self.pdf_toggle_btn.setStyleSheet(S.btn_ghost())
        self.pdf_toggle_btn.clicked.connect(self._toggle_pdf_settings)
        self.pdf_toggle_btn.hide()
        layout.addWidget(self.pdf_toggle_btn)

        self.settings_widget = QWidget()
        settings_panel = QVBoxLayout(self.settings_widget)
        settings_panel.setSpacing(10)
        settings_panel.setContentsMargins(20, 0, 20, 10)

        settings_row = QHBoxLayout()
        settings_row.setSpacing(15)

        # Resolution
        res_group = QGroupBox(self.tr("Resolution"))
        res_group.setFont(_font(S.SIZE_SMALL, bold=True))
        res_group.setStyleSheet(S.groupbox_style())
        res_layout = QVBoxLayout(res_group)
        self.resolution_combo = QComboBox()
        self.resolution_combo.setFont(_font())
        self.resolution_combo.setStyleSheet(S.combo_style())
        self.resolution_combo.addItems([
            self.tr("Coarse (~0.2 A)"),
            self.tr("Medium (~0.1 A)"),
            self.tr("Fine (~0.05 A)"),
        ])
        self.resolution_combo.setCurrentIndex(1)
        res_layout.addWidget(self.resolution_combo)
        settings_row.addWidget(res_group)

        # Slices
        slice_group = QGroupBox(self.tr("Slices"))
        slice_group.setFont(_font(S.SIZE_SMALL, bold=True))
        slice_group.setStyleSheet(S.groupbox_style())
        slice_layout = QHBoxLayout(slice_group)
        num_label = QLabel(self.tr("Number:"))
        num_label.setFont(_font())
        num_label.setStyleSheet(S.label_title())
        slice_layout.addWidget(num_label)
        self.slice_spinbox = QSpinBox()
        self.slice_spinbox.setFont(_font())
        self.slice_spinbox.setStyleSheet(S.spinbox_style())
        self.slice_spinbox.setRange(5, 30)
        self.slice_spinbox.setValue(15)
        slice_layout.addWidget(self.slice_spinbox)
        settings_row.addWidget(slice_group)

        # Output
        output_group = QGroupBox(self.tr("Output"))
        output_group.setFont(_font(S.SIZE_SMALL, bold=True))
        output_group.setStyleSheet(S.groupbox_style())
        output_layout = QVBoxLayout(output_group)
        self.promolecule_check = QCheckBox(self.tr("Promolecule"))
        self.promolecule_check.setFont(_font())
        self.promolecule_check.setStyleSheet(S.checkbox_style())
        self.promolecule_check.setChecked(True)
        output_layout.addWidget(self.promolecule_check)
        self.deformation_check = QCheckBox(self.tr("Deformation"))
        self.deformation_check.setFont(_font())
        self.deformation_check.setStyleSheet(S.checkbox_style())
        self.deformation_check.setChecked(True)
        output_layout.addWidget(self.deformation_check)
        settings_row.addWidget(output_group)

        # Color Mode
        color_group = QGroupBox(self.tr("Color Mode"))
        color_group.setFont(_font(S.SIZE_SMALL, bold=True))
        color_group.setStyleSheet(S.groupbox_style())
        color_layout = QVBoxLayout(color_group)
        self.color_combo = QComboBox()
        self.color_combo.setFont(_font())
        self.color_combo.setStyleSheet(S.combo_style())
        self.color_combo.addItems([
            self.tr("Black & White"),
            self.tr("Color"),
        ])
        color_layout.addWidget(self.color_combo)
        settings_row.addWidget(color_group)

        # Detail Level
        detail_group = QGroupBox(self.tr("Detail Level"))
        detail_group.setFont(_font(S.SIZE_SMALL, bold=True))
        detail_group.setStyleSheet(S.groupbox_style())
        detail_layout = QVBoxLayout(detail_group)
        self.detail_combo = QComboBox()
        self.detail_combo.setFont(_font())
        self.detail_combo.setStyleSheet(S.combo_style())
        self.detail_combo.addItems([
            self.tr("Simple"),
            self.tr("Advanced"),
        ])
        self.detail_combo.setCurrentIndex(0)
        detail_layout.addWidget(self.detail_combo)
        settings_row.addWidget(detail_group)

        settings_panel.addLayout(settings_row)

        # Generate button
        gen_row = QHBoxLayout()
        gen_row.addStretch()
        self.generate_btn = QPushButton(self.tr("Generate PDFs"))
        self.generate_btn.setFont(_font(13, bold=True))
        self.generate_btn.setStyleSheet(S.btn_primary())
        self.generate_btn.clicked.connect(self._on_generate)
        gen_row.addWidget(self.generate_btn)
        gen_row.addStretch()
        settings_panel.addLayout(gen_row)

        self.settings_widget.hide()
        layout.addWidget(self.settings_widget)

        return screen

    # =========================================================================
    # Processing Screen
    # =========================================================================

    def _create_processing_screen(self) -> QWidget:
        screen = QWidget()
        screen.setObjectName("processing_screen")
        screen.setStyleSheet("#processing_screen { background-color: #ffffff; }")
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setSpacing(20)

        layout.addStretch()

        title = QLabel(self.tr("Generating..."))
        title.setFont(_font(24, bold=True))
        title.setStyleSheet(S.label_title())
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setMinimumWidth(400)
        self.progress_bar.setStyleSheet(S.progress_bar())
        layout.addWidget(self.progress_bar, alignment=Qt.AlignmentFlag.AlignCenter)

        self.progress_label = QLabel(self.tr("Starting..."))
        self.progress_label.setFont(_font())
        self.progress_label.setStyleSheet(S.label_secondary())
        self.progress_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.progress_label)

        layout.addStretch()

        return screen

    # =========================================================================
    # Complete Screen
    # =========================================================================

    def _create_complete_screen(self) -> QWidget:
        screen = QWidget()
        screen.setObjectName("complete_screen")
        screen.setStyleSheet("#complete_screen { background-color: #ffffff; }")
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setSpacing(20)

        layout.addStretch()

        title = QLabel(self.tr("Complete!"))
        title.setFont(_font(28, bold=True))
        title.setStyleSheet(S.label_success())
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        self.complete_summary = QLabel()
        self.complete_summary.setFont(_font())
        self.complete_summary.setStyleSheet(S.label_secondary())
        self.complete_summary.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.complete_summary)

        layout.addSpacing(20)

        button_layout = QHBoxLayout()
        button_layout.addStretch()

        self.open_folder_btn = QPushButton(self.tr("Open Folder"))
        self.open_folder_btn.setFont(_font())
        self.open_folder_btn.setStyleSheet(S.btn_secondary())
        self.open_folder_btn.clicked.connect(self._open_output_folder)
        button_layout.addWidget(self.open_folder_btn)

        self.new_molecule_btn = QPushButton(self.tr("New Molecule"))
        self.new_molecule_btn.setFont(_font())
        self.new_molecule_btn.setStyleSheet(S.btn_primary())
        self.new_molecule_btn.clicked.connect(self._go_home)
        button_layout.addWidget(self.new_molecule_btn)

        button_layout.addStretch()
        layout.addLayout(button_layout)

        layout.addStretch()

        return screen

    # =========================================================================
    # Logic (unchanged)
    # =========================================================================

    def _cleanup_worker(self, worker: QThread | None) -> None:
        if worker is not None and worker.isRunning():
            worker.quit()
            worker.wait(1000)
            if worker.isRunning():
                worker.terminate()

    def _on_search(self):
        search_text = self.search_input.text().strip()
        if not search_text:
            return

        self.search_button.setEnabled(False)
        self.search_input.setEnabled(False)
        self.status_label.setText(self.tr("Searching..."))

        self._cleanup_worker(self.search_worker)

        self.search_worker = SearchWorker(self.searcher, search_text)
        self.search_worker.finished.connect(self._on_search_complete)
        self.search_worker.start()

    def _on_search_complete(self, results: list[MoleculeResult], is_online: bool):
        self.search_button.setEnabled(True)
        self.search_input.setEnabled(True)
        self.status_label.setText("")

        self.current_results = results
        query = self.search_input.text().strip()

        self.results_title.setText(self.tr(f'Results for "{query}"'))

        if is_online:
            self.online_status.setText(self.tr("COD Online"))
            self.online_status.setStyleSheet(f"color: {S.SUCCESS};")
        else:
            self.online_status.setText(self.tr("Offline (local only)"))
            self.online_status.setStyleSheet(f"color: {S.ACCENT};")

        while self.results_layout.count() > 1:
            child = self.results_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        if results:
            for mol in results:
                card = self._create_molecule_card(mol)
                self.results_layout.insertWidget(self.results_layout.count() - 1, card)
        else:
            no_results = QLabel(self.tr("No molecules found"))
            no_results.setFont(_font())
            no_results.setAlignment(Qt.AlignmentFlag.AlignCenter)
            no_results.setStyleSheet(f"color: {S.TEXT_MUTED}; padding: 40px;")
            self.results_layout.insertWidget(0, no_results)

        self.stacked_widget.setCurrentWidget(self.results_screen)

    def _create_molecule_card(self, mol: MoleculeResult) -> QWidget:
        card = QFrame()
        card.setFrameShape(QFrame.Shape.StyledPanel)
        card.setStyleSheet(S.card_style())

        layout = QHBoxLayout(card)

        info = QVBoxLayout()

        name = QLabel(mol.name)
        name.setFont(_font(S.SIZE_BODY, bold=True))
        name.setStyleSheet(S.label_title())
        info.addWidget(name)

        source = "Local" if mol.is_local else f"COD {mol.cod_id}"
        details = QLabel(f"{mol.formula} · {mol.atom_count} atoms · {source}")
        details.setFont(_font())
        details.setStyleSheet(S.label_secondary())
        info.addWidget(details)

        if mol.space_group:
            sg = QLabel(f"Space group: {mol.space_group}")
            sg.setFont(_font(S.SIZE_SMALL))
            sg.setStyleSheet(S.label_muted())
            info.addWidget(sg)

        layout.addLayout(info)
        layout.addStretch()

        select_btn = QPushButton(self.tr("Select"))
        select_btn.setFont(_font())
        select_btn.setStyleSheet(S.btn_secondary())
        select_btn.clicked.connect(lambda: self._on_select_molecule(mol))
        layout.addWidget(select_btn)

        return card

    def _on_select_molecule(self, mol: MoleculeResult):
        self.selected_molecule = mol
        self.selected_structure = self.searcher.get_structure(mol)

        self.preview_title.setText(mol.name)
        self.mol_info_label.setText(
            f"{mol.formula} · {mol.atom_count} atoms · Space group: {mol.space_group}"
        )

        if self.selected_structure:
            self.molecule_viewer.set_structure(self.selected_structure)
            self._setup_selection_manager()
        else:
            self.molecule_viewer.clear()

        self.slice_explorer.clear()
        self._switch_view("2d")
        self.stacked_widget.setCurrentWidget(self.preview_screen)

    def _setup_selection_manager(self):
        if not self.selected_structure:
            return

        from dorothy.core.selection import SelectionManager
        import numpy as np

        coords = self.selected_structure.get_cartesian_coords(
            align_to_principal_axes=True
        )

        self.selection_manager = SelectionManager(self)
        self.selection_manager.set_coordinates(coords)
        self.selection_manager.selection_changed.connect(self._on_selection_changed)
        self.selection_manager.plane_defined.connect(self._on_plane_defined)

        self.molecule_viewer.canvas.atom_picked.connect(
            self.selection_manager.toggle_atom
        )
        self.slice_explorer.canvas.atom_picked.connect(
            self.selection_manager.toggle_atom
        )

        self.molecule_viewer.canvas.set_pick_enabled(True)
        self.slice_explorer.canvas.set_pick_enabled(True)

    def _on_selection_changed(self, indices: list):
        self.molecule_viewer.canvas.set_selection(indices)
        self.slice_explorer.canvas.set_selection(indices)

        if len(indices) == 0:
            hint = "Click 4 atoms: 3 for plane + 1 for orientation"
        elif len(indices) < 3:
            hint = f"Selected {len(indices)}/4 atoms (need 3 for plane)"
        elif len(indices) == 3:
            hint = f"Selected 3/4 atoms (select 1 more for 'up' direction)"
        else:
            hint = "Plane defined! Slices centered on selected atoms."

        if self.selected_molecule:
            self.mol_info_label.setText(
                f"{self.selected_molecule.formula} - {hint}"
            )

    def _on_plane_defined(self, plane_definition):
        self._current_plane_definition = plane_definition
        self.slice_explorer._promolecule_cube = None
        self.slice_explorer._deformation_cube = None
        self.slice_explorer.canvas._density_cube = None
        if self.selected_structure:
            self._start_3d_density_calculation()

    def _clear_atom_selection(self):
        if self.selection_manager:
            self.selection_manager.clear()

    def _reset_slice_plane(self):
        self._clear_atom_selection()
        self._current_plane_definition = None
        self.slice_explorer._promolecule_cube = None
        self.slice_explorer._deformation_cube = None
        self.slice_explorer.canvas._density_cube = None
        if self.selected_structure:
            self._start_3d_density_calculation()

    def _on_generate(self):
        if not self.selected_structure:
            QMessageBox.warning(self, "Error", "No molecule structure available.")
            return

        wants_deformation = self.deformation_check.isChecked()
        if wants_deformation and not is_xtb_installed():
            dialog = XtbDownloadDialog(self)
            result = dialog.exec()
            if result == QDialog.DialogCode.Rejected:
                self.deformation_check.setChecked(False)
            elif not is_xtb_installed():
                self.deformation_check.setChecked(False)

        default_name = self.selected_molecule.name.lower().replace(' ', '_')
        output_dir = QFileDialog.getExistingDirectory(
            self,
            self.tr("Select Output Directory"),
            str(Path.home() / "Desktop"),
        )

        if not output_dir:
            return

        output_path = Path(output_dir) / f"{default_name}_dorothy"

        resolution_idx = self.resolution_combo.currentIndex()
        resolution = ["coarse", "medium", "fine"][resolution_idx]
        detail_level = "simple" if self.detail_combo.currentIndex() == 0 else "advanced"

        settings = GenerationSettings(
            resolution=resolution,
            n_slices=self.slice_spinbox.value(),
            generate_promolecule=self.promolecule_check.isChecked(),
            generate_deformation=self.deformation_check.isChecked(),
            color_mode="bw" if self.color_combo.currentIndex() == 0 else "color",
            detail_level=detail_level,
        )

        self.progress_bar.setValue(0)
        self.progress_label.setText(self.tr("Starting..."))
        self.stacked_widget.setCurrentWidget(self.processing_screen)

        self._cleanup_worker(self.generation_worker)

        self.generation_worker = GenerationWorker(
            self.selected_structure, output_path, settings
        )
        self.generation_worker.progress.connect(self._on_generation_progress)
        self.generation_worker.finished.connect(self._on_generation_complete)
        self.generation_worker.start()

    def _on_generation_progress(self, message: str, percent: int):
        self.progress_bar.setValue(percent)
        self.progress_label.setText(message)

    def _on_generation_complete(self, result: GenerationResult):
        self.last_result = result

        if result.success:
            n_promol = len(result.promolecule_pdfs)
            n_deform = len(result.deformation_pdfs)
            total = n_promol + n_deform

            summary = f"Generated {total} PDF slices\n"
            if n_promol > 0:
                summary += f"Promolecule: {n_promol} slices\n"
            if n_deform > 0:
                summary += f"Deformation: {n_deform} slices\n"
            if result.used_xtb:
                summary += "\n(Used xTB for molecular density)"
            else:
                summary += "\n(Promolecule density only - install xTB for deformation)"

            self.complete_summary.setText(summary)
            self.stacked_widget.setCurrentWidget(self.complete_screen)
        else:
            QMessageBox.critical(
                self,
                "Generation Failed",
                f"Error: {result.error_message}"
            )
            self.stacked_widget.setCurrentWidget(self.preview_screen)

    def _open_output_folder(self):
        if self.last_result and self.last_result.output_dir:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.last_result.output_dir)))

    def _go_home(self):
        self.stacked_widget.setCurrentWidget(self.home_screen)

    def _go_to_results(self):
        self.stacked_widget.setCurrentWidget(self.results_screen)

    def _toggle_pdf_settings(self):
        visible = self.settings_widget.isVisible()
        self.settings_widget.setVisible(not visible)
        if visible:
            self.pdf_toggle_btn.setText(self.tr("PDF Export Settings  >"))
        else:
            self.pdf_toggle_btn.setText(self.tr("PDF Export Settings  v"))

    def _switch_view(self, view: str):
        if view == "2d":
            self.view_2d_btn.setChecked(True)
            self.view_3d_btn.setChecked(False)
            self.viewer_stack.setCurrentWidget(self.molecule_viewer)
            self.pdf_toggle_btn.hide()
            self.settings_widget.hide()
        else:
            self.view_2d_btn.setChecked(False)
            self.view_3d_btn.setChecked(True)
            self.viewer_stack.setCurrentWidget(self.slice_explorer)
            self.pdf_toggle_btn.show()

            if self.selected_structure and self.slice_explorer._promolecule_cube is None:
                self._start_3d_density_calculation()

    def _start_3d_density_calculation(self):
        if not self.selected_structure:
            return

        self.mol_info_label.setText("Calculating density...")
        self._cleanup_worker(self.xtb_density_worker)

        plane_def = getattr(self, '_current_plane_definition', None)

        self.xtb_density_worker = XtbDensityWorker(
            self.selected_structure,
            plane_definition=plane_def
        )
        self.xtb_density_worker.progress.connect(self._on_xtb_density_progress)
        self.xtb_density_worker.finished.connect(self._on_xtb_density_complete)
        self.xtb_density_worker.start()

    def _on_xtb_density_progress(self, message: str):
        self.mol_info_label.setText(message)

    def _on_xtb_density_complete(self, promolecule, molecular, deformation):
        if promolecule:
            n_slices = self.slice_spinbox.value()
            self.slice_explorer.set_density_cubes(
                promolecule=promolecule,
                molecular=molecular,
                deformation=deformation,
                n_slices=n_slices
            )

            if deformation:
                self.mol_info_label.setText(
                    f"{self.selected_molecule.formula} - "
                    "Molecular & deformation density available"
                )
            else:
                self.mol_info_label.setText(
                    f"{self.selected_molecule.formula} - "
                    "Promolecule only (install xTB for molecular)"
                )

    def _update_3d_preview(self):
        if not self.selected_structure:
            return

        cube = create_density_cube_from_structure(
            self.selected_structure,
            resolution="coarse",
            align_to_principal_axes=True
        )

        n_slices = self.slice_spinbox.value()
        self.slice_explorer.set_density_cubes(
            promolecule=cube,
            deformation=None,
            n_slices=n_slices
        )

    def tr(self, text: str) -> str:
        return QCoreApplication.translate("MainWindow", text)
