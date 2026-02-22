# Dorothy: Crystallography Teaching Tool

**Honoring Dorothy Hodgkin's pioneering work in X-ray crystallography.**

Dorothy recreates historical electron density contour maps as a hands-on teaching tool. Students physically stack transparent slices to visualize 3D electron density and deduce molecular structure—essentially reverse-engineering what pioneering crystallographers did in the 1930s–1960s.

## Current Status

**Version:** 0.7.0

### What's Working

- **PyQt6 Desktop Application** - Native UI with Fusion style, design system
- **Molecule Search** - Local CIF files + COD online (with fallback) + load any CIF file
- **CIF Parser** - Extracts atomic coordinates, cell parameters, space groups, symmetry expansion
- **2D Molecule Viewer** - Bond detection, CPK coloring, interactive atom selection, 3D rotation
- **3D Slice Explorer** - Interactive 3D preview with zoom, bonds, and custom orientation
- **Interactive Plane Selection** - Click 4 atoms to define custom slicing plane (3 for plane + 1 for orientation)
- **Measurement Tool** - Click atoms to measure distances (2), angles (3), and dihedrals (4)
- **Auto xTB on 3D View** - Automatically calculates deformation density when switching to 3D
- **Generation Pipeline** - Promolecule density calculation, contour slicing, PDF export
- **Principal Axes Alignment** - Molecules auto-rotated for optimal slicing orientation
- **xTB Integration** - Deformation density via xTB (Homebrew/conda), with install dialog
- **Heatmap Mode** - Color gradient rendering for screen visualization (vs contours for print)
- **Detail Level Combo** - Fast (0.20 Å) / Good (0.10 Å) / Best (0.05 Å) integrated in toolbar
- **Planes Toggle** - Show/hide slice plane rectangles to view contours in isolation
- **Help Dialog** - Comprehensive "What is Dorothy?" modal with scientific background and references
- **Fixed Contour Levels** - Advanced mode uses consistent levels across slices for comparability
- **i18n Ready** - Qt translation system in place (English only for now)

### What's Not Yet Complete

- Settings persistence

---

## Quick Start

```bash
# Clone/download the project
cd Dorothy

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -e .

# Run
python -m dorothy.main
```

### Dependencies

- Python 3.10+
- PyQt6
- NumPy
- Matplotlib
- Requests

---

## Project Structure

```
Dorothy/
├── dorothy/
│   ├── __init__.py
│   ├── main.py                 # Application entry point
│   ├── core/
│   │   ├── constants.py        # Physical constants and settings
│   │   ├── cif_parser.py       # CIF file parsing
│   │   ├── cod_search.py       # COD + local molecule search
│   │   ├── density.py          # Electron density calculations
│   │   ├── contours.py         # Contour generation, PDF export
│   │   ├── generator.py        # Main generation pipeline
│   │   ├── selection.py        # Atom selection for plane definition
│   │   └── xtb_manager.py      # xTB download and execution
│   ├── ui/
│   │   ├── main_window.py      # Main application window
│   │   ├── molecule_viewer.py  # 2D structure viewer
│   │   ├── slice_explorer.py   # 3D interactive slice preview
│   │   └── styles.py           # Design system (colors, fonts, widget styles)
│   └── resources/
│       └── translations/       # i18n files
├── tests/                      # pytest test suite
│   ├── test_density.py         # Density calculation tests
│   ├── test_cif_parser.py      # CIF parsing tests
│   └── test_selection.py       # Plane selection tests
├── examples/
│   └── aspirin/                # Example CIF files
├── pyproject.toml
└── readme.md
```

---

## The Two Cubes Concept

1. **Promolecule Density** — Superposition of spherical, non-interacting atoms. Shows "here are the atoms."

2. **Deformation Density** — Molecular density minus promolecule density. Isolates bonding effects: accumulation in bonds/lone pairs, depletion elsewhere.

This pairing teaches both structure determination and chemical bonding concepts.

### Detail Levels

| Mode | Contour Levels | Best For |
|------|----------------|----------|
| **Simple** | Auto-scaled per slice | Beginners, cleaner output |
| **Advanced** | Fixed levels (0.01, 0.02, 0.04...) | Shows σ and π bonds |

### Contour Styling

| Density Type | Positive | Negative |
|--------------|----------|----------|
| **Promolecule** | Solid black lines | — |
| **Deformation** | Solid (blue) = accumulation | Dashed (red) = depletion |

---

## Adding Example Molecules

Place CIF files in `examples/<molecule_name>/`:

```
examples/
├── aspirin/
│   └── 185472.cif
├── benzene/
│   └── <any>.cif
└── caffeine/
    └── <any>.cif
```

The app auto-discovers these on startup.

---

## Next Steps

### Medium Priority
1. Add more bundled example molecules
2. Settings persistence (QSettings)
3. ZIP export option

### Future
4. App packaging for distribution (macOS .dmg, Windows .exe, Linux AppImage)
5. Additional languages (Portuguese, Spanish)

---

## Technical Notes

### Why xTB?

- Open source (LGPL-3.0), lightweight
- Handles organic molecules well
- Fast (seconds for small molecules)
- Auto-downloaded, not bundled (clean licensing)

### Deformation Density Calculation

The deformation density uses coordinates directly from the xTB cube file output to ensure perfect alignment between molecular and promolecule densities. By default, molecules are aligned to principal axes before xTB calculation so that Z-slices cut parallel to the molecular plane.

When using custom plane selection (4 atoms), the molecule is rotated so the user-selected plane becomes horizontal. This triggers a full density recalculation with the rotated coordinates, ensuring contours are computed directly on the selected plane orientation—not interpolated from a tilted slice.

### Grid Resolution

| Setting | Spacing | Use Case |
|---------|---------|----------|
| Coarse | ~0.2 Å | Quick preview |
| Medium | ~0.1 Å | Default teaching use |
| Fine | ~0.05 Å | Publication quality |

### Physical Output

- Default 15 slices per cube (adjustable 5-30)
- A5 size (148 × 210 mm)
- Print on transparent film, apply to PVC sheets, stack to visualize 3D

---

## Historical Context

Dorothy Hodgkin (1910–1994) pioneered X-ray crystallography, determining structures of penicillin, vitamin B12, and insulin. She drew electron density contours on transparent Perspex sheets and stacked them to visualize molecular structure.

### Hodgkin's workflow vs Dorothy's workflow

Hodgkin's original method ran in one direction:

> Crystal → X-ray diffraction → phase problem → electron density map → atom positions

She began with a physical crystal, collected a diffraction pattern, solved the phase problem, and used the resulting electron density to locate atoms she did not yet know.

Dorothy (the application) runs in the **opposite** direction:

> Atom positions (from CIF) → computed electron density → visualization

We already know where the atoms are — that information was determined by crystallographers and deposited in the COD. Starting from those known positions, Dorothy computes the electron density (either as a promolecule superposition or via xTB quantum chemistry) and visualizes it.

Neither density shown by Dorothy comes from an X-ray experiment directly. What comes from the experiment are the atom positions in the CIF file. The densities are computed from those positions, making them theoretical reconstructions of what a crystallographer would have observed in an experimental density map.

In a teaching context, this reversal is powerful: students see the electron density that Hodgkin would have painstakingly reconstructed from diffraction data, but without needing a crystal, an X-ray source, or months of hand calculation.

---

## Citation

If you use Dorothy in your work, please cite:

> https://doi.org/10.5281/zenodo.18733780

## License

GPL-3.0 — see [LICENSE](LICENSE)

xTB is LGPL-3.0 (downloaded separately, not bundled)

---

*February 2025*

---

## Changelog

### v0.7.0 (February 2025)
- **Symmetry Expansion**: CIF parser now applies space group symmetry operations to reconstruct complete molecules from asymmetric units, with automatic molecule extraction via bond connectivity
- **Disorder Resolution**: Partial-occupancy atoms are grouped by element and reduced to the correct count using occupancy sums (e.g. 8 F atoms at 0.25 occupancy → keep 2 representatives), correctly reconstructing disordered structures like As2F11
- **Extended Element Support**: Added As, Se, B, Si to element tables (colors, radii, density params)
- **Load CIF File**: New "load a CIF file" link on home screen to open any CIF file directly (e.g. from CCDC, ICSD, or custom sources)
- **Cross-platform Builds**: Standalone executables for macOS, Windows, and Linux via PyInstaller + GitHub Actions
- **xTB Download Fixes**: Corrected download URLs for Windows, removed non-existent macOS binaries
- **Citation**: Added Zenodo DOI citation to help dialog and readme

### v0.6.0 (February 2025)
- **WYSIWYG Viewer**: 3D viewer controls now match PDF export settings — what you see is what you get
- **Detail Level Combo**: Replaced standalone resolution spinbox with Fast/Good/Best dropdown integrated into the toolbar
- **Measurement Tool**: Click atoms to measure distances (2 atoms), angles (3 atoms), and dihedral angles (4 atoms) with gold overlay annotations
- **Planes Toggle**: New checkbox to show/hide slice plane rectangles while keeping contour lines visible
- **Help Dialog**: Comprehensive "What is Dorothy?" modal covering the science, workflow, Dorothy Hodgkin's legacy, and full academic references
- **Workflow Explanation**: Help dialog explains how Dorothy reverses Hodgkin's original crystallographic workflow (structure→density vs density→structure)
- **Design System**: Centralized styles module with Fusion rendering for consistent cross-platform appearance
- **Full grid slicing**: Z-slices now span the full density grid so density fades in from zero and back out
- **Flexible resolution**: `create_density_cube_from_structure` accepts both string presets and numeric Å values
- **Bond density labels**: Clearer bond annotation format (`density: 0.123` instead of `ρ=0.123`)

### v0.5.1 (December 2024)
- **Centralized constants**: New `constants.py` module for physical constants, element data, and settings
- **Test suite**: Added 42 pytest tests for density calculations, CIF parsing, and plane selection
- **Proper logging**: Replaced `print()` statements with Python `logging` module
- **Worker thread cleanup**: Fixed potential thread accumulation in background workers
- **Code quality**: Removed magic numbers, improved type hints

### v0.5.0 (December 2024)
- **Heatmap Rendering Mode**: Toggle between contour lines and color gradient display
  - Contour mode: Traditional line-based visualization (best for printing)
  - Heatmap mode: Color-mapped bitmap with `RdYlBu_r` colormap (best for screen)
- **2D Preview supports heatmaps**: View 2D button shows heatmap with colorbar when in heatmap mode

### v0.4.2 (December 2024)
- **Logo on home screen**: Replaced text title with Dorothy logo image
- **Improved COD reliability**: Added retry logic with exponential backoff for COD API searches and CIF downloads
- **Cleaner UI flow**: Export settings now only appear in 3D view mode, keeping 2D view focused on molecule visualization

### v0.4.1 (December 2024)
- **Fixed animation crash**: Promolecule and deformation cubes now use matching grids, preventing ValueError when blending densities during animation

### v0.4.0 (December 2024)
- **View 2D Button**: Opens current slice as a 2D contour plot with atom positions overlaid
- **Animate Button**: Smooth animation morphing from promolecule → molecular density
  - Shows how electron density changes when bonds form
  - Helps students understand deformation density concept
- **Click on Bonds**: Click any bond to see length and local density
  - Highlighted bonds show distance (Å) and density (e/Bohr³)
  - Click again to deselect
- **Improved Slice Range**: Slices now cover molecule extent (0.5Å padding) instead of full cube
- **Density at Point**: New `get_density_at_point()` API for querying density values
- **2D Preview shows atoms**: Atoms within 0.5Å of slice are overlaid on contour plot

### v0.3.8 (December 2024)
- **Improved 3D Explorer UI**: Reorganized into cleaner two-row layout
  - Row 1: View options (density type, color mode, display toggles, spacing)
  - Row 2: Navigation controls (slice slider, zoom, rotation, reset)
- **Z-coordinate display**: Slice label now shows position in Angstroms (e.g., "8/15 (z=1.2Å)")
- **Color mode toggle buttons**: Replaced dropdown with B&W/Color toggle buttons
- **Reset View button**: Single button resets both zoom and rotation to defaults
- **Tooltips everywhere**: All controls now have helpful tooltips
- Renamed confusing "Contours" slider to "Spacing"

### v0.3.7 (December 2024)
- **Single Slice Toggle**: New checkbox in 3D view to show only the selected slice, hiding all other planes and contours
- Easier to focus on one slice at a time without visual clutter

### v0.3.6 (December 2024)
- **2D Molecule Rotation**: Drag mouse or use arrow buttons to rotate molecule in 2D view
- Depth-based rendering gives pseudo-3D effect (size/alpha varies with depth)
- **Fixed COD CIF parsing**: Parser now correctly stops at loop boundaries, preventing atom data from being mixed with anisotropic displacement parameters and geometry data
- COD-downloaded molecules now display with correct coordinates and bond lengths

### v0.3.5 (December 2024)
- **Fixed 3D view distortion**: Added equal aspect ratio to prevent molecule stretching
- **Fixed custom plane coordinate mismatch**: Plane rotation now uses consistent coordinate system
- Custom plane selection properly aligns selected atoms to XY plane
- **Fixed COD search**: Updated API to use POST with text1 parameter, handle SSL certificate mismatch
- **Fixed COD CIF parsing**: Download and parse CIF files from COD, estimate atom count from formula
- Fixed space group parsing for quoted values in CIF files

### v0.3.4 (December 2024)
- **Architecture fix**: Custom plane selection now rotates the molecule before density calculation
- Plane selection triggers full xTB recalculation with rotated coordinates
- Contours are computed directly on the selected plane, not interpolated
- Removed complex oriented slice interpolation code (~200 lines)
- Fixed density info display: proper unit conversion (e/Bohr³ → e/Å³)
- Code cleanup: removed unused imports and dead code

### v0.3.3 (December 2024)
- Major fix: Custom plane slices now display as horizontal planes (like auto mode)
- Slices distributed over molecule extent, not full density cube
- Both auto and plane-selection modes now show consistent horizontal slices
- Oriented slices sample from tilted planes but display horizontally for comparison

### v0.3.2 (December 2024)
- Fixed contour rendering for custom plane selection (matplotlib 3.8+ compatibility)
- Added bilinear interpolation for smooth contour lines on oriented slices
- Fixed slider minimum width to prevent UI deformation
- Improved contour-to-3D coordinate mapping accuracy

### v0.3.1 (December 2024)
- 4-atom plane selection: 3 atoms define plane, 4th defines "up" direction, center is centroid
- Contour level slider: adjust contour spacing (tighter/looser)
- Density info display: toggle to show min/max/total electron density
- Improved density normalization using proper volume integration
- Fixed contour rendering with valid level filtering

### v0.3.0 (December 2024)
- Interactive plane selection: click atoms in 2D or 3D view to define custom slicing orientation
- Auto-run xTB when switching to 3D view (with progress indicator)
- Zoom controls for 3D cube view (+/- buttons, scroll wheel, reset)
- Rotation controls for 3D view (arrow buttons + drag-to-rotate)
- Bond visualization in 3D view (toggle on/off)
- Selection synchronization between 2D and 3D views
- Custom plane slicing via `get_oriented_slices()` with scipy interpolation

### v0.2.1 (December 2024)
- Added 3D Interactive Slice Explorer with contour visualization
- Toggle between 2D structure and 3D slice preview
- Slice navigation slider with highlighted current slice
- Mouse rotation/zoom for 3D view
- B&W and color mode support
