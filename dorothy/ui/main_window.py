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
from PyQt6.QtGui import QFont, QDesktopServices, QPixmap
from PyQt6.QtCore import QUrl

from dorothy.core.cod_search import CODSearch, MoleculeResult
from dorothy.core.cif_parser import MoleculeStructure
from dorothy.core.generator import GenerationPipeline, GenerationSettings, GenerationResult
from dorothy.core.xtb_manager import is_xtb_installed, download_xtb, get_download_url, get_xtb_install_instructions
from dorothy.core.density import create_density_cube_from_structure
from dorothy.ui.molecule_viewer import MoleculeViewer
from dorothy.ui.slice_explorer import SliceExplorer


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
    finished = pyqtSignal(object, object, object)  # promolecule_cube, deformation_cube, elf_cube

    def __init__(self, structure: MoleculeStructure, plane_definition=None, calculate_elf: bool = True):
        super().__init__()
        self.structure = structure
        self.plane_definition = plane_definition
        self.calculate_elf = calculate_elf

    def run(self):
        from dorothy.core.xtb_manager import is_xtb_installed, run_xtb_density
        from dorothy.core.density import (
            parse_cube_file,
            calculate_deformation_density,
            create_density_cube_from_structure
        )
        import tempfile

        promolecule = None
        deformation = None
        elf_cube = None

        try:
            # Always create promolecule (fast) as fallback
            # Use plane rotation if provided, otherwise use principal axes
            self.progress.emit("Generating promolecule density...")
            promolecule = create_density_cube_from_structure(
                self.structure,
                resolution="coarse",
                align_to_principal_axes=(self.plane_definition is None),
                plane_definition=self.plane_definition
            )

            # Try xTB if available
            if is_xtb_installed():
                self.progress.emit("Running xTB calculation...")
                with tempfile.TemporaryDirectory() as tmp_dir:
                    result = run_xtb_density(
                        self.structure,
                        Path(tmp_dir),
                        lambda msg: self.progress.emit(msg),
                        plane_definition=self.plane_definition,
                        calculate_elf=self.calculate_elf
                    )

                    # Handle both old (single path) and new (tuple) return formats
                    if self.calculate_elf:
                        cube_path, elf_path = result
                    else:
                        cube_path = result
                        elf_path = None

                    if cube_path:
                        self.progress.emit("Processing density data...")
                        molecular = parse_cube_file(cube_path)
                        if molecular:
                            # Get both promolecule and deformation on the same grid
                            # This ensures they can be blended for animation
                            promolecule, deformation = calculate_deformation_density(
                                molecular,
                                self.structure
                            )

                    if elf_path:
                        self.progress.emit("Processing ELF data...")
                        elf_cube = parse_cube_file(elf_path)

            self.progress.emit("Ready")

        except Exception as e:
            self.progress.emit(f"Error: {e}")

        self.finished.emit(promolecule, deformation, elf_cube)


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

        # Title
        title = QLabel("xTB Not Installed")
        title_font = QFont()
        title_font.setPointSize(14)
        title_font.setBold(True)
        title.setFont(title_font)
        layout.addWidget(title)

        # Explanation
        explanation = QLabel(
            "xTB is required for deformation density calculation.\n\n"
            "Without xTB, only promolecule density (atomic positions) "
            "will be generated. With xTB, you also get deformation density "
            "showing chemical bonding."
        )
        explanation.setWordWrap(True)
        explanation.setStyleSheet("color: #666;")
        layout.addWidget(explanation)

        # Installation instructions
        instructions = get_xtb_install_instructions()
        install_label = QLabel(instructions)
        install_label.setWordWrap(True)
        install_label.setStyleSheet(
            "background-color: #f5f5f5; padding: 10px; "
            "font-family: monospace; border-radius: 4px;"
        )
        install_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(install_label)

        # Progress bar (for auto-download, hidden initially)
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.hide()
        layout.addWidget(self.progress_bar)

        # Status label
        self.status_label = QLabel()
        self.status_label.setStyleSheet("color: #666;")
        self.status_label.hide()
        layout.addWidget(self.status_label)

        # Buttons
        button_layout = QHBoxLayout()

        self.skip_btn = QPushButton("Skip (Promolecule Only)")
        self.skip_btn.clicked.connect(self.reject)
        button_layout.addWidget(self.skip_btn)

        if self.can_auto_download:
            self.download_btn = QPushButton("Download xTB")
            self.download_btn.setStyleSheet("""
                QPushButton {
                    background-color: #2563eb;
                    color: white;
                    font-weight: bold;
                    padding: 8px 16px;
                    border-radius: 4px;
                }
                QPushButton:hover {
                    background-color: #1d4ed8;
                }
            """)
            self.download_btn.clicked.connect(self._start_download)
            button_layout.addWidget(self.download_btn)
        else:
            # No auto-download available, just show "I've Installed It" button
            self.check_btn = QPushButton("I've Installed It")
            self.check_btn.setStyleSheet("""
                QPushButton {
                    background-color: #2563eb;
                    color: white;
                    font-weight: bold;
                    padding: 8px 16px;
                    border-radius: 4px;
                }
                QPushButton:hover {
                    background-color: #1d4ed8;
                }
            """)
            self.check_btn.clicked.connect(self._check_installation)
            button_layout.addWidget(self.check_btn)

        layout.addLayout(button_layout)

    def _check_installation(self):
        """Check if user has installed xTB."""
        if is_xtb_installed():
            self.accept()
        else:
            self.status_label.show()
            self.status_label.setText("xTB not found. Please install it and try again.")
            self.status_label.setStyleSheet("color: #dc2626;")

    def _start_download(self):
        """Start the xTB download."""
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
        """Handle download progress."""
        if total > 0:
            percent = int(100 * downloaded / total)
            self.progress_bar.setValue(percent)
            mb_downloaded = downloaded / (1024 * 1024)
            mb_total = total / (1024 * 1024)
            self.status_label.setText(f"Downloading... {mb_downloaded:.1f} / {mb_total:.1f} MB")

    def _on_finished(self, success: bool):
        """Handle download completion."""
        if success:
            self.status_label.setText("Download complete!")
            self.status_label.setStyleSheet("color: #16a34a;")
            self.accept()
        else:
            self.status_label.setText("Download failed. Please install manually.")
            self.status_label.setStyleSheet("color: #dc2626;")
            self.download_btn.setEnabled(True)
            self.skip_btn.setEnabled(True)


class MainWindow(QMainWindow):
    """Main application window for Dorothy."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Dorothy")
        self.setMinimumSize(800, 600)

        self.searcher = CODSearch()  # Single instance to cache local molecules
        self.search_worker = None
        self.generation_worker = None
        self.xtb_density_worker = None
        self.current_results = []
        self.selected_molecule: MoleculeResult | None = None
        self.selected_structure: MoleculeStructure | None = None
        self.last_result: GenerationResult | None = None

        # Selection manager for plane definition
        self.selection_manager = None

        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Create stacked widget for different screens
        self.stacked_widget = QStackedWidget()
        main_layout.addWidget(self.stacked_widget)

        # Add screens
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

        # Show home screen
        self.stacked_widget.setCurrentWidget(self.home_screen)

    def _create_home_screen(self) -> QWidget:
        """Create the home/search screen."""
        screen = QWidget()
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setSpacing(20)

        # Logo
        logo_label = QLabel()
        logo_path = Path(__file__).parent.parent.parent / "logo" / "dorothy_logo.png"
        if logo_path.exists():
            pixmap = QPixmap(str(logo_path))
            # Scale to reasonable size while keeping aspect ratio
            scaled = pixmap.scaledToWidth(400, Qt.TransformationMode.SmoothTransformation)
            logo_label.setPixmap(scaled)
        else:
            # Fallback to text if logo not found
            logo_label.setText(self.tr("Dorothy"))
            title_font = QFont()
            title_font.setPointSize(32)
            title_font.setBold(True)
            logo_label.setFont(title_font)
        logo_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(logo_label)

        # Subtitle
        subtitle = QLabel(self.tr("Crystallography Teaching Tool"))
        subtitle_font = QFont()
        subtitle_font.setPointSize(12)
        subtitle.setFont(subtitle_font)
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        subtitle.setStyleSheet("color: #666;")
        layout.addWidget(subtitle)

        layout.addSpacing(40)

        # Search section
        search_container = QWidget()
        search_layout = QVBoxLayout(search_container)
        search_layout.setSpacing(10)

        # Search input
        search_input_layout = QHBoxLayout()
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText(self.tr("Search molecule..."))
        self.search_input.setMinimumHeight(40)
        search_input_font = QFont()
        search_input_font.setPointSize(12)
        self.search_input.setFont(search_input_font)
        self.search_input.returnPressed.connect(self._on_search)
        search_input_layout.addWidget(self.search_input)

        self.search_button = QPushButton(self.tr("Search"))
        self.search_button.setMinimumHeight(40)
        self.search_button.setMinimumWidth(100)
        self.search_button.clicked.connect(self._on_search)
        search_input_layout.addWidget(self.search_button)

        search_layout.addLayout(search_input_layout)

        # Status label (for showing searching... or errors)
        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.status_label.setStyleSheet("color: #666;")
        search_layout.addWidget(self.status_label)

        # Example molecules
        examples = QLabel(
            self.tr("Examples: aspirin · benzene · caffeine · urea · naphthalene")
        )
        examples.setAlignment(Qt.AlignmentFlag.AlignCenter)
        examples.setStyleSheet("color: #888;")
        search_layout.addWidget(examples)

        layout.addWidget(search_container)

        # Add stretch to push content to top
        layout.addStretch()

        # Footer
        footer = QLabel(
            self.tr(
                "Honoring Dorothy Hodgkin's pioneering work in X-ray crystallography"
            )
        )
        footer.setAlignment(Qt.AlignmentFlag.AlignCenter)
        footer.setStyleSheet("color: #999; font-size: 10pt;")
        layout.addWidget(footer)

        return screen

    def _create_results_screen(self) -> QWidget:
        """Create the search results screen."""
        screen = QWidget()
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)

        # Header with back button
        header = QHBoxLayout()
        back_btn = QPushButton(self.tr("< Back"))
        back_btn.clicked.connect(self._go_home)
        back_btn.setMaximumWidth(80)
        header.addWidget(back_btn)

        self.results_title = QLabel(self.tr("Results"))
        title_font = QFont()
        title_font.setPointSize(18)
        title_font.setBold(True)
        self.results_title.setFont(title_font)
        header.addWidget(self.results_title)
        header.addStretch()

        self.online_status = QLabel()
        header.addWidget(self.online_status)

        layout.addLayout(header)

        # Scrollable results area
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

    def _create_preview_screen(self) -> QWidget:
        """Create the molecule preview screen with settings."""
        screen = QWidget()
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(15)

        # Header with back button
        header = QHBoxLayout()
        back_btn = QPushButton(self.tr("< Back"))
        back_btn.clicked.connect(self._go_to_results)
        back_btn.setMaximumWidth(80)
        header.addWidget(back_btn)

        self.preview_title = QLabel(self.tr("Molecule"))
        title_font = QFont()
        title_font.setPointSize(18)
        title_font.setBold(True)
        self.preview_title.setFont(title_font)
        header.addWidget(self.preview_title)
        header.addStretch()

        layout.addLayout(header)

        # Main content: viewer + settings side by side
        content = QHBoxLayout()

        # Left: Molecule viewer with view toggle
        viewer_container = QVBoxLayout()

        # View toggle buttons
        view_toggle = QHBoxLayout()
        self.view_2d_btn = QPushButton(self.tr("2D Structure"))
        self.view_2d_btn.setCheckable(True)
        self.view_2d_btn.setChecked(True)
        self.view_2d_btn.clicked.connect(lambda: self._switch_view("2d"))
        self.view_2d_btn.setStyleSheet("""
            QPushButton {
                padding: 8px 16px;
                border: 1px solid #ccc;
                border-radius: 4px 0 0 4px;
                background: #2563eb;
                color: white;
                font-weight: bold;
            }
            QPushButton:!checked {
                background: #f0f0f0;
                color: #333;
            }
        """)
        view_toggle.addWidget(self.view_2d_btn)

        self.view_3d_btn = QPushButton(self.tr("3D Slices"))
        self.view_3d_btn.setCheckable(True)
        self.view_3d_btn.clicked.connect(lambda: self._switch_view("3d"))
        self.view_3d_btn.setStyleSheet("""
            QPushButton {
                padding: 8px 16px;
                border: 1px solid #ccc;
                border-left: none;
                border-radius: 0 4px 4px 0;
                background: #f0f0f0;
                color: #333;
            }
            QPushButton:checked {
                background: #2563eb;
                color: white;
                font-weight: bold;
            }
        """)
        view_toggle.addWidget(self.view_3d_btn)

        # Selection controls for plane definition
        view_toggle.addSpacing(20)
        self.clear_selection_btn = QPushButton(self.tr("Clear Selection"))
        self.clear_selection_btn.setMaximumWidth(120)
        self.clear_selection_btn.clicked.connect(self._clear_atom_selection)
        self.clear_selection_btn.setToolTip("Clear atom selection for plane definition")
        view_toggle.addWidget(self.clear_selection_btn)

        self.reset_plane_btn = QPushButton(self.tr("Reset Plane"))
        self.reset_plane_btn.setMaximumWidth(100)
        self.reset_plane_btn.clicked.connect(self._reset_slice_plane)
        self.reset_plane_btn.setToolTip("Reset to default slicing orientation")
        view_toggle.addWidget(self.reset_plane_btn)

        view_toggle.addStretch()
        viewer_container.addLayout(view_toggle)

        # Stacked widget for 2D/3D views
        self.viewer_stack = QStackedWidget()

        # 2D Molecule viewer
        self.molecule_viewer = MoleculeViewer()
        self.molecule_viewer.setMinimumSize(350, 350)
        self.viewer_stack.addWidget(self.molecule_viewer)

        # 3D Slice explorer
        self.slice_explorer = SliceExplorer()
        self.slice_explorer.setMinimumSize(350, 350)
        self.viewer_stack.addWidget(self.slice_explorer)

        viewer_container.addWidget(self.viewer_stack)

        # Molecule info below viewer
        self.mol_info_label = QLabel()
        self.mol_info_label.setStyleSheet("color: #666;")
        self.mol_info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        viewer_container.addWidget(self.mol_info_label)

        content.addLayout(viewer_container, stretch=1)

        # Right: Settings panel (only visible in 3D view)
        self.settings_widget = QWidget()
        settings_panel = QVBoxLayout(self.settings_widget)
        settings_panel.setSpacing(15)
        settings_panel.setContentsMargins(0, 0, 0, 0)

        # Resolution setting
        res_group = QGroupBox(self.tr("Resolution"))
        res_layout = QVBoxLayout(res_group)
        self.resolution_combo = QComboBox()
        self.resolution_combo.addItems([
            self.tr("Coarse (~0.2 A) - Quick preview"),
            self.tr("Medium (~0.1 A) - Default"),
            self.tr("Fine (~0.05 A) - High quality"),
        ])
        self.resolution_combo.setCurrentIndex(1)
        res_layout.addWidget(self.resolution_combo)
        settings_panel.addWidget(res_group)

        # Slice settings
        slice_group = QGroupBox(self.tr("Slices"))
        slice_layout = QHBoxLayout(slice_group)
        slice_layout.addWidget(QLabel(self.tr("Number of slices:")))
        self.slice_spinbox = QSpinBox()
        self.slice_spinbox.setRange(5, 30)
        self.slice_spinbox.setValue(15)
        slice_layout.addWidget(self.slice_spinbox)
        settings_panel.addWidget(slice_group)

        # Output options
        output_group = QGroupBox(self.tr("Output"))
        output_layout = QVBoxLayout(output_group)
        self.promolecule_check = QCheckBox(self.tr("Promolecule density"))
        self.promolecule_check.setChecked(True)
        output_layout.addWidget(self.promolecule_check)
        self.deformation_check = QCheckBox(self.tr("Deformation density"))
        self.deformation_check.setChecked(True)
        output_layout.addWidget(self.deformation_check)
        settings_panel.addWidget(output_group)

        # Color mode
        color_group = QGroupBox(self.tr("Color Mode"))
        color_layout = QVBoxLayout(color_group)
        self.color_combo = QComboBox()
        self.color_combo.addItems([
            self.tr("Black & White"),
            self.tr("Color (Blue/Red)"),
        ])
        color_layout.addWidget(self.color_combo)
        settings_panel.addWidget(color_group)

        # Detail level
        detail_group = QGroupBox(self.tr("Detail Level"))
        detail_layout = QVBoxLayout(detail_group)
        self.detail_combo = QComboBox()
        self.detail_combo.addItems([
            self.tr("Simple (for beginners)"),
            self.tr("Advanced (shows π-bonds)"),
        ])
        self.detail_combo.setCurrentIndex(0)
        detail_layout.addWidget(self.detail_combo)
        detail_hint = QLabel(self.tr("Advanced mode uses fixed contour levels\nto reveal subtle bonding features."))
        detail_hint.setStyleSheet("color: #888; font-size: 9pt;")
        detail_hint.setWordWrap(True)
        detail_layout.addWidget(detail_hint)
        settings_panel.addWidget(detail_group)

        settings_panel.addStretch()

        # Generate button
        self.generate_btn = QPushButton(self.tr("Generate PDFs"))
        self.generate_btn.setMinimumHeight(50)
        self.generate_btn.setStyleSheet("""
            QPushButton {
                background-color: #2563eb;
                color: white;
                font-size: 14pt;
                font-weight: bold;
                border-radius: 8px;
            }
            QPushButton:hover {
                background-color: #1d4ed8;
            }
            QPushButton:pressed {
                background-color: #1e40af;
            }
        """)
        self.generate_btn.clicked.connect(self._on_generate)
        settings_panel.addWidget(self.generate_btn)

        # Start hidden (shown when switching to 3D view)
        self.settings_widget.hide()
        content.addWidget(self.settings_widget)

        layout.addLayout(content)

        return screen

    def _create_processing_screen(self) -> QWidget:
        """Create the processing/progress screen."""
        screen = QWidget()
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setSpacing(20)

        layout.addStretch()

        # Title
        title = QLabel(self.tr("Generating..."))
        title_font = QFont()
        title_font.setPointSize(24)
        title_font.setBold(True)
        title.setFont(title_font)
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setMinimumWidth(400)
        layout.addWidget(self.progress_bar, alignment=Qt.AlignmentFlag.AlignCenter)

        # Status message
        self.progress_label = QLabel(self.tr("Starting..."))
        self.progress_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.progress_label.setStyleSheet("color: #666;")
        layout.addWidget(self.progress_label)

        layout.addStretch()

        return screen

    def _create_complete_screen(self) -> QWidget:
        """Create the completion screen."""
        screen = QWidget()
        layout = QVBoxLayout(screen)
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setSpacing(20)

        layout.addStretch()

        # Success icon/title
        title = QLabel(self.tr("Complete!"))
        title_font = QFont()
        title_font.setPointSize(28)
        title_font.setBold(True)
        title.setFont(title_font)
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("color: #16a34a;")
        layout.addWidget(title)

        # Summary
        self.complete_summary = QLabel()
        self.complete_summary.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.complete_summary.setStyleSheet("font-size: 12pt; color: #666;")
        layout.addWidget(self.complete_summary)

        layout.addSpacing(20)

        # Buttons
        button_layout = QHBoxLayout()
        button_layout.addStretch()

        self.open_folder_btn = QPushButton(self.tr("Open Folder"))
        self.open_folder_btn.setMinimumHeight(40)
        self.open_folder_btn.setMinimumWidth(120)
        self.open_folder_btn.clicked.connect(self._open_output_folder)
        button_layout.addWidget(self.open_folder_btn)

        self.new_molecule_btn = QPushButton(self.tr("New Molecule"))
        self.new_molecule_btn.setMinimumHeight(40)
        self.new_molecule_btn.setMinimumWidth(120)
        self.new_molecule_btn.clicked.connect(self._go_home)
        button_layout.addWidget(self.new_molecule_btn)

        button_layout.addStretch()
        layout.addLayout(button_layout)

        layout.addStretch()

        return screen

    def _on_search(self):
        """Handle search button click."""
        search_text = self.search_input.text().strip()
        if not search_text:
            return

        # Disable search while running
        self.search_button.setEnabled(False)
        self.search_input.setEnabled(False)
        self.status_label.setText(self.tr("Searching..."))

        # Run search in background thread
        self.search_worker = SearchWorker(self.searcher, search_text)
        self.search_worker.finished.connect(self._on_search_complete)
        self.search_worker.start()

    def _on_search_complete(self, results: list[MoleculeResult], is_online: bool):
        """Handle search completion."""
        self.search_button.setEnabled(True)
        self.search_input.setEnabled(True)
        self.status_label.setText("")

        self.current_results = results
        query = self.search_input.text().strip()

        # Update results title
        self.results_title.setText(self.tr(f'Results for "{query}"'))

        # Update online status
        if is_online:
            self.online_status.setText(self.tr("COD Online"))
            self.online_status.setStyleSheet("color: #2a2; font-size: 10pt;")
        else:
            self.online_status.setText(self.tr("Offline (local only)"))
            self.online_status.setStyleSheet("color: #a52; font-size: 10pt;")

        # Clear previous results
        while self.results_layout.count() > 1:
            child = self.results_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        # Add new results
        if results:
            for mol in results:
                card = self._create_molecule_card(mol)
                self.results_layout.insertWidget(self.results_layout.count() - 1, card)
        else:
            no_results = QLabel(self.tr("No molecules found"))
            no_results.setAlignment(Qt.AlignmentFlag.AlignCenter)
            no_results.setStyleSheet("color: #888; padding: 40px;")
            self.results_layout.insertWidget(0, no_results)

        # Show results screen
        self.stacked_widget.setCurrentWidget(self.results_screen)

    def _create_molecule_card(self, mol: MoleculeResult) -> QWidget:
        """Create a card widget for a molecule result."""
        card = QFrame()
        card.setFrameShape(QFrame.Shape.StyledPanel)
        card.setStyleSheet("""
            QFrame {
                background: #f8f8f8;
                border: 1px solid #ddd;
                border-radius: 8px;
                padding: 10px;
            }
            QFrame:hover {
                background: #f0f0f0;
                border-color: #bbb;
            }
        """)

        layout = QHBoxLayout(card)

        # Info section
        info = QVBoxLayout()

        name = QLabel(mol.name)
        name_font = QFont()
        name_font.setPointSize(12)
        name_font.setBold(True)
        name.setFont(name_font)
        info.addWidget(name)

        # Source badge (Local or COD ID)
        source = "Local" if mol.is_local else f"COD {mol.cod_id}"
        details = QLabel(f"{mol.formula} · {mol.atom_count} atoms · {source}")
        details.setStyleSheet("color: #666;")
        info.addWidget(details)

        if mol.space_group:
            sg = QLabel(f"Space group: {mol.space_group}")
            sg.setStyleSheet("color: #888; font-size: 10pt;")
            info.addWidget(sg)

        layout.addLayout(info)
        layout.addStretch()

        # Select button
        select_btn = QPushButton(self.tr("Select"))
        select_btn.setMinimumWidth(80)
        select_btn.clicked.connect(lambda: self._on_select_molecule(mol))
        layout.addWidget(select_btn)

        return card

    def _on_select_molecule(self, mol: MoleculeResult):
        """Handle molecule selection - show preview screen."""
        self.selected_molecule = mol

        # Get structure data
        self.selected_structure = self.searcher.get_structure(mol)

        # Update preview screen
        self.preview_title.setText(mol.name)
        self.mol_info_label.setText(
            f"{mol.formula} · {mol.atom_count} atoms · Space group: {mol.space_group}"
        )

        # Show structure in viewer
        if self.selected_structure:
            self.molecule_viewer.set_structure(self.selected_structure)

            # Set up selection manager for plane definition
            self._setup_selection_manager()
        else:
            self.molecule_viewer.clear()

        # Clear 3D preview (will be regenerated when switching to 3D view)
        self.slice_explorer.clear()

        # Reset to 2D view
        self._switch_view("2d")

        # Show preview screen
        self.stacked_widget.setCurrentWidget(self.preview_screen)

    def _setup_selection_manager(self):
        """Set up selection manager for the current structure."""
        if not self.selected_structure:
            return

        from dorothy.core.selection import SelectionManager
        import numpy as np

        # Get aligned coordinates for plane calculation
        coords = self.selected_structure.get_cartesian_coords(
            align_to_principal_axes=True
        )

        # Create selection manager
        self.selection_manager = SelectionManager(self)
        self.selection_manager.set_coordinates(coords)

        # Connect selection signals
        self.selection_manager.selection_changed.connect(self._on_selection_changed)
        self.selection_manager.plane_defined.connect(self._on_plane_defined)

        # Connect picker signals from both views
        self.molecule_viewer.canvas.atom_picked.connect(
            self.selection_manager.toggle_atom
        )
        self.slice_explorer.canvas.atom_picked.connect(
            self.selection_manager.toggle_atom
        )

        # Enable picking by default (can be toggled)
        self.molecule_viewer.canvas.set_pick_enabled(True)
        self.slice_explorer.canvas.set_pick_enabled(True)

    def _on_selection_changed(self, indices: list):
        """Synchronize selection across views."""
        self.molecule_viewer.canvas.set_selection(indices)
        self.slice_explorer.canvas.set_selection(indices)

        # Update status with 4-atom selection info
        if len(indices) == 0:
            hint = "Click 4 atoms: 3 for plane + 1 for orientation"
        elif len(indices) < 3:
            hint = f"Selected {len(indices)}/4 atoms (need 3 for plane)"
        elif len(indices) == 3:
            hint = f"Selected 3/4 atoms (select 1 more for 'up' direction)"
        else:
            hint = "Plane defined! Slices centered on selected atoms."

        # Show hint in info label temporarily
        if self.selected_molecule:
            self.mol_info_label.setText(
                f"{self.selected_molecule.formula} - {hint}"
            )

    def _on_plane_defined(self, plane_definition):
        """Handle plane definition from 4 selected atoms.

        When a new plane is defined, we need to recalculate the density cube
        with the molecule rotated so the selected plane is horizontal.
        This ensures Z-slices cut parallel to the user-selected plane.
        """
        # Store the plane definition
        self._current_plane_definition = plane_definition

        # Clear cached densities to force recalculation
        self.slice_explorer._promolecule_cube = None
        self.slice_explorer._deformation_cube = None
        self.slice_explorer.canvas._density_cube = None

        # Trigger recalculation with the new plane orientation
        if self.selected_structure:
            self._start_3d_density_calculation()

    def _clear_atom_selection(self):
        """Clear current atom selection."""
        if self.selection_manager:
            self.selection_manager.clear()

    def _reset_slice_plane(self):
        """Reset to default principal-axes slicing."""
        self._clear_atom_selection()
        self._current_plane_definition = None

        # Clear cached densities to force recalculation with default orientation
        self.slice_explorer._promolecule_cube = None
        self.slice_explorer._deformation_cube = None
        self.slice_explorer.canvas._density_cube = None

        # Trigger recalculation with default (principal axes) orientation
        if self.selected_structure:
            self._start_3d_density_calculation()

    def _on_generate(self):
        """Handle generate PDFs button click."""
        if not self.selected_structure:
            QMessageBox.warning(self, "Error", "No molecule structure available.")
            return

        # Check if xTB is needed and not installed
        wants_deformation = self.deformation_check.isChecked()
        if wants_deformation and not is_xtb_installed():
            # Show install dialog
            dialog = XtbDownloadDialog(self)
            result = dialog.exec()

            if result == QDialog.DialogCode.Rejected:
                # User chose to skip - disable deformation
                self.deformation_check.setChecked(False)
            elif not is_xtb_installed():
                # Installation not complete
                self.deformation_check.setChecked(False)

        # Ask for output directory
        default_name = self.selected_molecule.name.lower().replace(' ', '_')
        output_dir = QFileDialog.getExistingDirectory(
            self,
            self.tr("Select Output Directory"),
            str(Path.home() / "Desktop"),
        )

        if not output_dir:
            return

        output_path = Path(output_dir) / f"{default_name}_dorothy"

        # Get settings
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

        # Show processing screen
        self.progress_bar.setValue(0)
        self.progress_label.setText(self.tr("Starting..."))
        self.stacked_widget.setCurrentWidget(self.processing_screen)

        # Run generation in background
        self.generation_worker = GenerationWorker(
            self.selected_structure, output_path, settings
        )
        self.generation_worker.progress.connect(self._on_generation_progress)
        self.generation_worker.finished.connect(self._on_generation_complete)
        self.generation_worker.start()

    def _on_generation_progress(self, message: str, percent: int):
        """Handle generation progress updates."""
        self.progress_bar.setValue(percent)
        self.progress_label.setText(message)

    def _on_generation_complete(self, result: GenerationResult):
        """Handle generation completion."""
        self.last_result = result

        if result.success:
            # Update completion screen
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
        """Open the output folder in file browser."""
        if self.last_result and self.last_result.output_dir:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.last_result.output_dir)))

    def _go_home(self):
        """Return to home screen."""
        self.stacked_widget.setCurrentWidget(self.home_screen)

    def _go_to_results(self):
        """Return to results screen."""
        self.stacked_widget.setCurrentWidget(self.results_screen)

    def _switch_view(self, view: str):
        """Switch between 2D and 3D views."""
        if view == "2d":
            self.view_2d_btn.setChecked(True)
            self.view_3d_btn.setChecked(False)
            self.viewer_stack.setCurrentWidget(self.molecule_viewer)
            self.settings_widget.hide()
        else:
            self.view_2d_btn.setChecked(False)
            self.view_3d_btn.setChecked(True)
            self.viewer_stack.setCurrentWidget(self.slice_explorer)
            self.settings_widget.show()

            # Auto-run xTB calculation for 3D preview if not already done
            if self.selected_structure and self.slice_explorer._promolecule_cube is None:
                self._start_3d_density_calculation()

    def _start_3d_density_calculation(self):
        """Start background xTB calculation for 3D preview."""
        if not self.selected_structure:
            return

        # Show progress indicator
        self.mol_info_label.setText("Calculating density...")

        # Get current plane definition if any
        plane_def = getattr(self, '_current_plane_definition', None)

        # Start background calculation with plane rotation
        self.xtb_density_worker = XtbDensityWorker(
            self.selected_structure,
            plane_definition=plane_def
        )
        self.xtb_density_worker.progress.connect(self._on_xtb_density_progress)
        self.xtb_density_worker.finished.connect(self._on_xtb_density_complete)
        self.xtb_density_worker.start()

    def _on_xtb_density_progress(self, message: str):
        """Handle xTB density progress updates."""
        self.mol_info_label.setText(message)

    def _on_xtb_density_complete(self, promolecule, deformation, elf_cube=None):
        """Handle xTB density calculation completion."""
        if promolecule:
            n_slices = self.slice_spinbox.value()
            self.slice_explorer.set_density_cubes(
                promolecule=promolecule,
                deformation=deformation,
                elf=elf_cube,
                n_slices=n_slices
            )

            # Update info label
            if deformation and elf_cube:
                self.mol_info_label.setText(
                    f"{self.selected_molecule.formula} - "
                    "Deformation density + ELF available"
                )
            elif deformation:
                self.mol_info_label.setText(
                    f"{self.selected_molecule.formula} - "
                    "Deformation density available (showing bonds)"
                )
            else:
                self.mol_info_label.setText(
                    f"{self.selected_molecule.formula} - "
                    "Promolecule only (install xTB for deformation)"
                )

    def _update_3d_preview(self):
        """Generate density cube for 3D preview (legacy sync method)."""
        if not self.selected_structure:
            return

        # Use coarse resolution for quick preview
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
        """Translate text using Qt's translation system."""
        return QCoreApplication.translate("MainWindow", text)
