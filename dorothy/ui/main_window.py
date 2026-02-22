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
# Worker threads
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

    def __init__(self, structure: MoleculeStructure, plane_definition=None,
                 resolution: str | float = 0.2):
        super().__init__()
        self.structure = structure
        self.plane_definition = plane_definition
        self.resolution = resolution

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
                resolution=self.resolution,
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


class PdfExportWorker(QThread):
    """Background thread for PDF export from viewer state."""
    progress = pyqtSignal(str, int)  # message, percent
    finished = pyqtSignal(bool, str, str)  # success, error_message, output_dir

    def __init__(self, promolecule_cube, deformation_cube, output_dir: Path,
                 n_slices: int, color_mode: str, molecule_name: str, formula: str):
        super().__init__()
        self.promolecule_cube = promolecule_cube
        self.deformation_cube = deformation_cube
        self.output_dir = Path(output_dir)
        self.n_slices = n_slices
        self.color_mode = color_mode
        self.molecule_name = molecule_name
        self.formula = formula

    def run(self):
        from dorothy.core.contours import (
            ContourSettings,
            generate_promolecule_contours,
            generate_deformation_contours,
            create_info_page,
        )

        try:
            self.output_dir.mkdir(parents=True, exist_ok=True)

            settings = ContourSettings(
                n_slices=self.n_slices,
                color_mode=self.color_mode,
            )

            if self.promolecule_cube:
                self.progress.emit("Generating promolecule slices...", 20)
                promol_dir = self.output_dir / "promolecule"
                generate_promolecule_contours(
                    self.promolecule_cube, promol_dir, settings,
                    lambda curr, total, msg: self.progress.emit(msg, 20 + int(30 * curr / total))
                )

            if self.deformation_cube:
                self.progress.emit("Generating deformation slices...", 55)
                deform_dir = self.output_dir / "deformation"
                generate_deformation_contours(
                    self.deformation_cube, deform_dir, settings,
                    lambda curr, total, msg: self.progress.emit(msg, 55 + int(30 * curr / total))
                )

            self.progress.emit("Creating info page...", 90)
            create_info_page(
                self.molecule_name, self.formula, settings,
                self.output_dir / "info.pdf"
            )

            self.progress.emit("Export complete!", 100)
            self.finished.emit(True, "", str(self.output_dir))

        except Exception as e:
            self.finished.emit(False, str(e), "")


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
# Help Dialog
# =============================================================================

class HelpDialog(QDialog):
    """Comprehensive help dialog explaining Dorothy's purpose, workflow, and science."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(self.tr("Dorothy — Help"))
        self.setMinimumSize(750, 600)
        self.setModal(True)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 16)
        layout.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        content = QLabel()
        content.setWordWrap(True)
        content.setTextFormat(Qt.TextFormat.RichText)
        content.setOpenExternalLinks(True)
        content.setFont(_font(S.SIZE_SUBTITLE))
        content.setStyleSheet(f"color: {S.TEXT}; padding: 32px 40px;")
        content.setText(self._build_html())

        scroll.setWidget(content)
        layout.addWidget(scroll)

        close_btn = QPushButton(self.tr("Close"))
        close_btn.setFont(_font(S.SIZE_BODY))
        close_btn.setStyleSheet(S.btn_secondary())
        close_btn.clicked.connect(self.accept)
        layout.addWidget(close_btn, alignment=Qt.AlignmentFlag.AlignCenter)

    def _build_html(self) -> str:
        a = S.ACCENT
        s = S.TEXT_SECONDARY
        f = S.FONT_FAMILY

        return self.tr(
            '<div style="font-family: \'{f}\'; line-height: 1.7; font-size: 14pt;">'

            # --- What is Dorothy? ---
            '<h2 style="color: {a}; margin-top: 0;">What is Dorothy?</h2>'
            '<p>Dorothy is a teaching tool for visualizing electron density '
            'inside molecules. It takes a molecular structure from a '
            'crystallographic database, computes the electron density on a '
            'three-dimensional grid, and slices that grid into a stack of '
            'two-dimensional contour maps suitable for printing on '
            'transparency sheets. When stacked, these sheets form a tangible, '
            'physical model of the electron cloud.</p>'
            '<p>Electron density is the central quantity in structural '
            'chemistry. It describes the probability of finding an electron '
            'at any given point in space, and its shape governs how atoms '
            'bond, how molecules react, and what properties a material will '
            'have. Despite its importance, electron density remains abstract '
            'for most students because it cannot be seen directly. Dorothy '
            'makes it concrete.</p>'

            # --- How to use Dorothy ---
            '<h2 style="color: {a};">How to use Dorothy</h2>'
            '<p><b>1. Search.</b> Enter a molecule name on the home screen '
            '(e.g. aspirin, benzene, caffeine). Dorothy queries the '
            'Crystallography Open Database for matching crystal structures.</p>'
            '<p><b>2. Select.</b> Choose a result from the list. The '
            'molecular formula, atom count, and space group are displayed. '
            'Dorothy downloads the crystal structure file (CIF) and extracts '
            'a single molecule from the unit cell.</p>'
            '<p><b>3. Configure.</b> Set the grid resolution, the number of '
            'slices, the color scheme, and which density types to generate.</p>'
            '<p><b>4. Generate.</b> Dorothy computes the electron density and '
            'produces the PDF slices. Depending on molecule size and '
            'resolution, this takes a few seconds to a few minutes.</p>'
            '<p><b>5. Explore.</b> Use the interactive 3D viewer to browse '
            'slices, rotate the molecule, switch density types, and toggle '
            'between contour lines, bonds, and heatmap overlays.</p>'

            # --- Where does the data come from? ---
            '<h2 style="color: {a};">Where does the data come from?</h2>'
            '<p>Molecular structures are drawn from the Crystallography Open '
            'Database (COD), a free, open-access repository containing over '
            '500,000 experimentally determined crystal structures '
            '(<a href="https://doi.org/10.1093/nar/gkr900">'
            'Gra\u017eulis <i>et al.</i>, 2012</a>).</p>'
            '<p>Each structure was determined by X-ray diffraction. In this '
            'technique, X-rays are passed through a crystal, producing a '
            'diffraction pattern that encodes the positions of the atoms in '
            'the lattice. From this pattern, crystallographers reconstruct '
            'the three-dimensional arrangement of every atom. The result is '
            'stored in a Crystallographic Information File (CIF), which '
            'records the fractional coordinates of each atom, the unit cell '
            'parameters (the size and shape of the repeating unit), and the '
            'space group (the symmetry operations relating molecules within '
            'the crystal) '
            '(<a href="https://doi.org/10.1107/S010876739101067X">'
            'Hall, Allen &amp; Brown, 1991</a>).</p>'
            '<p>Dorothy reads the CIF, applies the symmetry operations to '
            'generate all atoms in the unit cell, and extracts a single '
            'molecule by identifying bonded clusters.</p>'

            # --- What is electron density? ---
            '<h2 style="color: {a};">What is electron density?</h2>'
            '<p>Electron density is a scalar field that assigns a value, '
            'measured in electrons per cubic angstrom (e/\u00c5\u00b3), to '
            'every point in space. Near an atomic nucleus the density is '
            'high, reflecting the concentration of core and valence '
            'electrons. Far from any atom it approaches zero. Between bonded '
            'atoms a ridge of elevated density connects the two nuclei, '
            'making the chemical bond directly visible in the topology of '
            'the field '
            '(<a href="https://global.oup.com/academic/product/'
            'atoms-in-molecules-9780198558651">Bader, 1990</a>).</p>'
            '<p>Contour maps are the standard representation of electron '
            'density. Each contour line connects points of equal density, '
            'analogous to altitude contours on a topographic map. Tightly '
            'spaced contours indicate a steep gradient; widely spaced '
            'contours indicate a plateau.</p>'

            # --- Three types of density ---
            '<h2 style="color: {a};">Three types of density</h2>'

            '<h3 style="color: {a};">Promolecule density</h3>'
            '<p>The promolecule model places each atom at its '
            'crystallographic position but treats every atom as an isolated, '
            'non-interacting sphere. The total density is simply the '
            'superposition of these free-atom densities. This model serves '
            'as a chemically meaningful baseline: it represents what the '
            'electron distribution would look like in the absence of any '
            'bonding interaction '
            '(<a href="https://doi.org/10.1021/j100401a010">'
            'Spackman &amp; Maslen, 1986</a>).</p>'

            '<h3 style="color: {a};">Molecular density</h3>'
            '<p>The molecular density is computed from quantum mechanics. '
            'Dorothy uses the GFN2-xTB method (Extended Tight-Binding), a '
            'semi-empirical approach developed by the Grimme group at the '
            'University of Bonn '
            '(<a href="https://doi.org/10.1021/acs.jctc.8b01176">'
            'Bannwarth, Ehlert &amp; Grimme, 2019</a>; '
            '<a href="https://doi.org/10.1021/acs.jctc.7b00118">'
            'Grimme, Bannwarth &amp; Shushkov, 2017</a>). '
            'xTB solves an approximate Schr\u00f6dinger equation to '
            'determine how electrons redistribute when atoms form bonds. '
            'The result captures bonding density, lone pairs, and '
            'delocalization effects that the promolecule model misses '
            'entirely.</p>'
            '<p>As a semi-empirical method, xTB is considerably faster than '
            'full density functional theory while remaining accurate enough '
            'for pedagogical visualization. For a comprehensive review of '
            'the method family, see '
            '<a href="https://doi.org/10.1002/wcms.1493">'
            'Bannwarth <i>et al.</i>, 2021</a>.</p>'
            '<p style="color: {s};"><i>Note: xTB must be installed '
            'separately. If it is not available, only the promolecule '
            'density is generated.</i></p>'

            '<h3 style="color: {a};">Deformation density</h3>'
            '<p>The deformation density is the difference between molecular '
            'and promolecule densities. By subtracting the non-interacting '
            'atomic background, it isolates the redistribution of electrons '
            'caused by chemical bonding '
            '(<a href="https://global.oup.com/academic/product/'
            'x-ray-charge-densities-and-chemical-bonding-9780195098235">'
            'Coppens, 1997</a>).</p>'
            '<p>Positive regions indicate accumulation of electron density. '
            'These appear between bonded atoms (shared electrons) and around '
            'lone pairs. Negative regions indicate depletion, where density '
            'has migrated to participate in bonding elsewhere.</p>'
            '<p>The deformation density is the most direct visualization of '
            'how and where atoms share electrons. It is, in effect, a '
            'fingerprint of chemical bonding.</p>'

            # --- What are the PDF slices? ---
            '<h2 style="color: {a};">What are the PDF slices?</h2>'
            '<p>The three-dimensional electron density is a continuous field '
            'filling the space around the molecule. To make it physically '
            'accessible, Dorothy sections it into parallel two-dimensional '
            'slices along the Z-axis.</p>'
            '<p>Each slice is exported as an A5-format PDF displaying '
            'contour lines of constant electron density, with atom positions '
            'and bond connections overlaid. The intended workflow is to print '
            'each slice on a transparency sheet, then stack all sheets with '
            'uniform spacing. Looking through the assembled stack reveals a '
            'three-dimensional picture of the electron density that students '
            'can hold, rotate, and examine from any angle.</p>'

            # --- The 3D viewer ---
            '<h2 style="color: {a};">The interactive 3D viewer</h2>'
            '<p>After generating the slices, Dorothy opens an interactive '
            'viewer for detailed exploration:</p>'
            '<ul style="margin-left: 20px;">'
            '<li><b>Slice slider.</b> Scroll through slices from bottom to '
            'top to observe how the density varies along the molecular '
            'axis.</li>'
            '<li><b>Density type.</b> Switch between promolecule, molecular, '
            'and deformation views.</li>'
            '<li><b>Contour / Heatmap.</b> Toggle between contour lines '
            '(suited for printing) and color heatmaps (suited for on-screen '
            'analysis).</li>'
            '<li><b>B&amp;W / Color.</b> Select monochrome or colored '
            'contours. In deformation mode, blue marks positive regions and '
            'red marks negative regions.</li>'
            '<li><b>Spacing.</b> Adjust the interval between contour levels. '
            'Tighter spacing resolves finer features; looser spacing produces '
            'a cleaner image.</li>'
            '<li><b>Bonds / Contours / Planes / Single / Info.</b> Toggle '
            'individual display elements.</li>'
            '<li><b>Zoom and Rotate.</b> Use the on-screen controls or drag '
            'directly on the 3D canvas.</li>'
            '<li><b>Animate.</b> Watch a smooth interpolation from '
            'promolecule to molecular density, making bonding effects visible '
            'in real time.</li>'
            '</ul>'

            # --- Dorothy Hodgkin ---
            '<h2 style="color: {a};">Who was Dorothy Hodgkin?</h2>'
            '<p>Dorothy Crowfoot Hodgkin (1910\u20131994) was a British '
            'chemist who pioneered the use of X-ray crystallography to '
            'determine the three-dimensional structures of biologically '
            'important molecules. Her most celebrated achievements were the '
            'structure determinations of penicillin (1945) and vitamin '
            'B<sub>12</sub> (1956), both of which were considered '
            'intractable problems at the time. She also contributed to early '
            'crystallographic studies of cholesterol and other sterols, and '
            'spent 35 years pursuing the structure of insulin, which she '
            'finally solved in 1969.</p>'
            '<p>She received the <b>Nobel Prize in Chemistry in 1964</b>, '
            'the third woman to do so after Marie Curie (1911) and '
            'Ir\u00e8ne Joliot-Curie (1935).</p>'
            '<p>Her electron density maps, laboriously computed by hand and '
            'drawn on stacked glass sheets, are the direct inspiration for '
            'this application. Dorothy brings that method into the digital '
            'age so that every student can experience the same insight: '
            'seeing atoms and bonds made visible.</p>'

            # --- Hodgkin's method vs Dorothy ---
            '<h2 style="color: {a};">How Dorothy relates to Hodgkin\u2019s method</h2>'
            '<p>Hodgkin\u2019s original workflow ran in one direction:</p>'
            '<p style="text-align: center; font-size: 13pt; color: {s};">'
            'Crystal \u2192 X-ray diffraction \u2192 phase problem '
            '\u2192 electron density map \u2192 atom positions</p>'
            '<p>She began with a physical crystal, collected a diffraction '
            'pattern, solved the notoriously difficult '
            '<a href="https://doi.org/10.1107/97809553602060000001">'
            'phase problem</a>, and used the resulting electron density to '
            'locate atoms she did not yet know.</p>'
            '<p>Dorothy (the application) runs in the <b>opposite</b> '
            'direction:</p>'
            '<p style="text-align: center; font-size: 13pt; color: {s};">'
            'Atom positions (from CIF) \u2192 computed electron density '
            '\u2192 visualization</p>'
            '<p>We already know where the atoms are \u2014 that information '
            'was determined by crystallographers and deposited in the '
            'Crystallography Open Database. Starting from those known '
            'positions, Dorothy computes the electron density '
            'computationally (either as a promolecule superposition or via '
            'xTB quantum chemistry) and visualizes it.</p>'
            '<p>Neither density shown by Dorothy comes from an X-ray '
            'experiment directly. What comes from the experiment are the '
            'atom positions in the CIF file. The densities are computed '
            'from those positions, making them theoretical reconstructions '
            'of what a crystallographer would have observed in an '
            'experimental density map. For molecules with accurate crystal '
            'structures, the computed and experimental densities agree '
            'closely \u2014 the bonding features, lone pairs, and '
            'deformation patterns are qualitatively the same.</p>'
            '<p>In a teaching context, this reversal is powerful: students '
            'see the electron density that Hodgkin would have painstakingly '
            'reconstructed from diffraction data, but without needing a '
            'crystal, an X-ray source, or months of hand calculation. They '
            'can focus on understanding what the density reveals about '
            'chemical bonding.</p>'

            # --- Data sources and methods ---
            '<h2 style="color: {a};">Data sources and methods</h2>'
            '<ul style="margin-left: 20px;">'
            '<li><b>Molecular structures:</b> Crystallography Open Database '
            '(COD)</li>'
            '<li><b>Quantum chemistry engine:</b> xTB by the Grimme '
            'group</li>'
            '<li><b>Atomic density data:</b> Thakkar and Koga tables of '
            'spherical atomic densities '
            '(<a href="https://doi.org/10.1007/BF02341696">'
            'Koga, 1997</a>)</li>'
            '</ul>'

            # --- References ---
            '<h2 style="color: {a};">References</h2>'
            '<p style="font-size: 12pt; line-height: 1.6;">'

            'Gra\u017eulis, S. <i>et al.</i> (2012). Crystallography Open '
            'Database (COD): an open-access collection of crystal structures '
            'and platform for world-wide collaboration. <i>Nucleic Acids '
            'Research</i>, 40(D1), D420\u2013D427. '
            '<a href="https://doi.org/10.1093/nar/gkr900">'
            'DOI: 10.1093/nar/gkr900</a><br><br>'

            'Hall, S.R., Allen, F.H. &amp; Brown, I.D. (1991). The '
            'crystallographic information file (CIF): a new standard archive '
            'file for crystallography. <i>Acta Crystallographica Section A</i>, '
            '47(6), 655\u2013685. '
            '<a href="https://doi.org/10.1107/S010876739101067X">'
            'DOI: 10.1107/S010876739101067X</a><br><br>'

            'Bader, R.F.W. (1990). <i>Atoms in Molecules: A Quantum Theory</i>. '
            'Oxford University Press. ISBN: 978-0-19-855865-1<br><br>'

            'Spackman, M.A. &amp; Maslen, E.N. (1986). Chemical properties '
            'from the promolecule. <i>Journal of Physical Chemistry</i>, '
            '90(10), 2020\u20132027. '
            '<a href="https://doi.org/10.1021/j100401a010">'
            'DOI: 10.1021/j100401a010</a><br><br>'

            'Coppens, P. (1997). <i>X-Ray Charge Densities and Chemical '
            'Bonding</i>. Oxford University Press. '
            'ISBN: 978-0-19-509823-5<br><br>'

            'Grimme, S., Bannwarth, C. &amp; Shushkov, P. (2017). A robust '
            'and accurate tight-binding quantum chemical method for structures, '
            'vibrational frequencies, and noncovalent interactions of large '
            'molecular systems parametrized for all spd-block elements '
            '(Z = 1\u201386). <i>Journal of Chemical Theory and Computation</i>, '
            '13(5), 1989\u20132009. '
            '<a href="https://doi.org/10.1021/acs.jctc.7b00118">'
            'DOI: 10.1021/acs.jctc.7b00118</a><br><br>'

            'Bannwarth, C., Ehlert, S. &amp; Grimme, S. (2019). GFN2-xTB: '
            'an accurate and broadly parametrized self-consistent tight-binding '
            'quantum chemical method with multipole electrostatics and '
            'density-dependent dispersion contributions. <i>Journal of '
            'Chemical Theory and Computation</i>, 15(3), 1652\u20131671. '
            '<a href="https://doi.org/10.1021/acs.jctc.8b01176">'
            'DOI: 10.1021/acs.jctc.8b01176</a><br><br>'

            'Bannwarth, C. <i>et al.</i> (2021). Extended tight-binding '
            'quantum chemistry methods. <i>WIREs Computational Molecular '
            'Science</i>, 11, e01493. '
            '<a href="https://doi.org/10.1002/wcms.1493">'
            'DOI: 10.1002/wcms.1493</a><br><br>'

            'Koga, T. (1997). Analytical Hartree-Fock electron densities '
            'for atoms He through Lr. <i>Theoretical Chemistry Accounts</i>, '
            '95, 113\u2013130. '
            '<a href="https://doi.org/10.1007/BF02341696">'
            'DOI: 10.1007/BF02341696</a><br><br>'

            'Koga, T., Kanayama, K., Watanabe, S. &amp; Thakkar, A.J. (1999). '
            'Analytical Hartree-Fock wave functions subject to cusp and '
            'asymptotic constraints: He to Xe, Li+ to Cs+, H\u2212 to I\u2212. '
            '<i>International Journal of Quantum Chemistry</i>, 71(6), '
            '491\u2013497. '
            '<a href="https://doi.org/10.1002/(SICI)1097-461X(1999)71:6'
            '%3C491::AID-QUA6%3E3.0.CO;2-T">'
            'DOI: 10.1002/(SICI)1097-461X(1999)71:6</a><br><br>'

            'Koga, T., Kanayama, K., Watanabe, T., Imai, T. &amp; '
            'Thakkar, A.J. (2000). Analytical Hartree-Fock wave functions '
            'for the atoms Cs to Lr. <i>Theoretical Chemistry Accounts</i>, '
            '104, 411\u2013413. '
            '<a href="https://doi.org/10.1007/s002140000150">'
            'DOI: 10.1007/s002140000150</a>'

            '</p>'
            '</div>'
        ).format(f=f, a=a, s=s)


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

        self.searcher = CODSearch()
        self.search_worker: SearchWorker | None = None
        self.xtb_density_worker: XtbDensityWorker | None = None
        self.pdf_export_worker: PdfExportWorker | None = None
        self.current_results = []
        self.selected_molecule: MoleculeResult | None = None
        self.selected_structure: MoleculeStructure | None = None
        self.selection_manager = None
        self._density_show_processing = False

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
        self.stacked_widget.addWidget(self.home_screen)
        self.stacked_widget.addWidget(self.results_screen)
        self.stacked_widget.addWidget(self.preview_screen)
        self.stacked_widget.addWidget(self.processing_screen)

        self.stacked_widget.setCurrentWidget(self.home_screen)

        # Floating "?" help button — always visible, top-right
        self.help_btn = QPushButton("?")
        self.help_btn.setParent(central_widget)
        self.help_btn.setObjectName("help_btn")
        self.help_btn.setFont(_font(S.SIZE_SUBTITLE, bold=True))
        self.help_btn.setFixedSize(36, 36)
        self.help_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.help_btn.setToolTip(self.tr("About Dorothy"))
        self.help_btn.setStyleSheet(f"""
            #help_btn {{
                background-color: white;
                color: {S.TEXT_SECONDARY};
                border: 2px solid {S.BORDER};
                border-radius: 18px;
                font-size: {S.SIZE_SUBTITLE}pt;
                font-family: "{S.FONT_FAMILY}";
                font-weight: bold;
            }}
            #help_btn:hover {{
                border-color: {S.ACCENT};
                color: {S.ACCENT};
            }}
        """)
        self.help_btn.raise_()
        self.help_btn.clicked.connect(self._show_help)

        # Connect viewer export signal
        self.slice_explorer.export_requested.connect(self._on_export_pdf)

        # Show maximized AFTER help button exists, so resizeEvent positions it
        self.showMaximized()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, 'help_btn'):
            self.help_btn.move(self.centralWidget().width() - 52, 12)

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
    # Preview Screen — WYSIWYG: Resolution + Generate, viewer has all controls
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

        # Detail level + Generate button
        self.detail_combo = QComboBox()
        self.detail_combo.setFont(_font())
        self.detail_combo.setStyleSheet(S.combo_style())
        self.detail_combo.addItems(["Fast", "Good", "Best"])
        self.detail_combo.setToolTip(
            "Fast = quick preview (0.20 \u00c5)\n"
            "Good = teaching quality (0.10 \u00c5)\n"
            "Best = high detail (0.05 \u00c5)"
        )
        view_toggle.addWidget(self.detail_combo)

        view_toggle.addSpacing(8)

        self.generate_btn = QPushButton(self.tr("Generate"))
        self.generate_btn.setFont(_font(bold=True))
        self.generate_btn.setStyleSheet(S.btn_primary())
        self.generate_btn.setToolTip("Compute electron density at the selected detail level")
        self.generate_btn.clicked.connect(self._on_generate)
        view_toggle.addWidget(self.generate_btn)

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
    # Logic
    # =========================================================================

    def _cleanup_worker(self, worker: QThread | None) -> None:
        if worker is not None and worker.isRunning():
            worker.quit()
            worker.wait(1000)
            if worker.isRunning():
                worker.terminate()

    def _show_help(self):
        dialog = HelpDialog(self)
        dialog.exec()

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

    # =========================================================================
    # Generate — computes density cubes and shows in 3D viewer
    # =========================================================================

    def _on_generate(self):
        """Compute density at the selected detail level and show in viewer."""
        if not self.selected_structure:
            return

        detail_map = {"Fast": 0.20, "Good": 0.10, "Best": 0.05}
        resolution = detail_map.get(self.detail_combo.currentText(), 0.20)
        self._start_3d_density_calculation(
            resolution=resolution,
            show_processing=True
        )

    def _start_3d_density_calculation(self, resolution: str | float = "coarse",
                                       show_processing: bool = False):
        """Start density computation in background.

        Args:
            resolution: Grid spacing (float in Å) or string preset
            show_processing: If True, show the processing screen during computation
        """
        if not self.selected_structure:
            return

        self._density_show_processing = show_processing

        if show_processing:
            self.progress_bar.setValue(0)
            self.progress_label.setText(self.tr("Computing density..."))
            self.stacked_widget.setCurrentWidget(self.processing_screen)
            self.generate_btn.setEnabled(False)
        else:
            self.mol_info_label.setText("Calculating density...")

        self._cleanup_worker(self.xtb_density_worker)

        plane_def = getattr(self, '_current_plane_definition', None)

        self.xtb_density_worker = XtbDensityWorker(
            self.selected_structure,
            plane_definition=plane_def,
            resolution=resolution
        )
        self.xtb_density_worker.progress.connect(self._on_xtb_density_progress)
        self.xtb_density_worker.finished.connect(self._on_xtb_density_complete)
        self.xtb_density_worker.start()

    def _on_xtb_density_progress(self, message: str):
        if self._density_show_processing:
            self.progress_label.setText(message)
            # Estimate progress from message content
            msg_lower = message.lower()
            if "promolecule" in msg_lower:
                self.progress_bar.setValue(20)
            elif "running xtb" in msg_lower or "xtb" in msg_lower:
                self.progress_bar.setValue(40)
            elif "processing" in msg_lower:
                self.progress_bar.setValue(70)
            elif "ready" in msg_lower:
                self.progress_bar.setValue(95)
        else:
            self.mol_info_label.setText(message)

    def _on_xtb_density_complete(self, promolecule, molecular, deformation):
        self.generate_btn.setEnabled(True)

        if promolecule:
            n_slices = self.slice_explorer.get_n_slices()
            self.slice_explorer.set_density_cubes(
                promolecule=promolecule,
                molecular=molecular,
                deformation=deformation,
                n_slices=n_slices
            )

        # If we were showing the processing screen, return to preview in 3D mode
        if self._density_show_processing:
            self._density_show_processing = False
            self.stacked_widget.setCurrentWidget(self.preview_screen)
            self._switch_view("3d")

        # Update status message
        if self.selected_molecule:
            if deformation:
                self.mol_info_label.setText(
                    f"{self.selected_molecule.formula} - "
                    "Molecular & deformation density available"
                )
            elif promolecule:
                self.mol_info_label.setText(
                    f"{self.selected_molecule.formula} - "
                    "Promolecule only (install xTB for molecular)"
                )

    # =========================================================================
    # Export PDF — generates PDFs from current viewer state
    # =========================================================================

    def _on_export_pdf(self):
        """Export PDFs using current viewer settings (WYSIWYG)."""
        promolecule = self.slice_explorer._promolecule_cube
        deformation = self.slice_explorer._deformation_cube

        if not promolecule and not deformation:
            QMessageBox.warning(
                self, self.tr("No Data"),
                self.tr("Generate density first before exporting.")
            )
            return

        # Ask for output directory
        default_name = self.selected_molecule.name.lower().replace(' ', '_') if self.selected_molecule else "molecule"
        output_dir = QFileDialog.getExistingDirectory(
            self,
            self.tr("Select Output Directory"),
            str(Path.home() / "Desktop"),
        )
        if not output_dir:
            return

        output_path = Path(output_dir) / f"{default_name}_dorothy"

        # Get viewer settings
        n_slices = self.slice_explorer.get_n_slices()
        color_mode = self.slice_explorer.canvas._color_mode

        # Show processing screen
        self.progress_bar.setValue(0)
        self.progress_label.setText(self.tr("Exporting PDFs..."))
        self.stacked_widget.setCurrentWidget(self.processing_screen)

        # Start export worker
        self._cleanup_worker(self.pdf_export_worker)
        self.pdf_export_worker = PdfExportWorker(
            promolecule_cube=promolecule,
            deformation_cube=deformation,
            output_dir=output_path,
            n_slices=n_slices,
            color_mode=color_mode,
            molecule_name=self.selected_molecule.name if self.selected_molecule else "Unknown",
            formula=self.selected_molecule.formula if self.selected_molecule else "",
        )
        self.pdf_export_worker.progress.connect(self._on_pdf_export_progress)
        self.pdf_export_worker.finished.connect(self._on_pdf_export_complete)
        self.pdf_export_worker.start()

    def _on_pdf_export_progress(self, message: str, percent: int):
        self.progress_label.setText(message)
        self.progress_bar.setValue(percent)

    def _on_pdf_export_complete(self, success: bool, error_message: str, output_dir: str):
        # Return to preview screen in 3D view
        self.stacked_widget.setCurrentWidget(self.preview_screen)
        self._switch_view("3d")

        if success:
            reply = QMessageBox.information(
                self, self.tr("Export Complete"),
                self.tr(f"PDFs exported to:\n{output_dir}\n\nOpen folder?"),
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes
            )
            if reply == QMessageBox.StandardButton.Yes:
                QDesktopServices.openUrl(QUrl.fromLocalFile(output_dir))
        else:
            QMessageBox.critical(
                self, self.tr("Export Failed"),
                self.tr(f"Error: {error_message}")
            )

    # =========================================================================
    # Navigation
    # =========================================================================

    def _go_home(self):
        self.stacked_widget.setCurrentWidget(self.home_screen)

    def _go_to_results(self):
        self.stacked_widget.setCurrentWidget(self.results_screen)

    def _switch_view(self, view: str):
        if view == "2d":
            self.view_2d_btn.setChecked(True)
            self.view_3d_btn.setChecked(False)
            self.viewer_stack.setCurrentWidget(self.molecule_viewer)
        else:
            self.view_2d_btn.setChecked(False)
            self.view_3d_btn.setChecked(True)
            self.viewer_stack.setCurrentWidget(self.slice_explorer)

            if self.selected_structure and self.slice_explorer._promolecule_cube is None:
                self._start_3d_density_calculation()

    def tr(self, text: str) -> str:
        return QCoreApplication.translate("MainWindow", text)
