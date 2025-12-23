# Dorothy: Crystallography Teaching Tool

**Honoring Dorothy Hodgkin's pioneering work in X-ray crystallography.**

Dorothy recreates historical electron density contour maps as a hands-on teaching tool. Students physically stack transparent slices to visualize 3D electron density and deduce molecular structure—essentially reverse-engineering what pioneering crystallographers did in the 1930s–1960s.

## Current Status

**Version:** 0.3.3 (Development)

### What's Working

- **PyQt6 Desktop Application** - Native UI with 5 screens
- **Molecule Search** - Local CIF files + COD online (with fallback)
- **CIF Parser** - Extracts atomic coordinates, cell parameters, space groups
- **2D Molecule Viewer** - Bond detection, CPK coloring, interactive atom selection
- **3D Slice Explorer** - Interactive 3D preview with zoom, bonds, and custom orientation
- **Interactive Plane Selection** - Click 4 atoms to define custom slicing plane (3 for plane + 1 for orientation)
- **Auto xTB on 3D View** - Automatically calculates deformation density when switching to 3D
- **Generation Pipeline** - Promolecule density calculation, contour slicing, PDF export
- **Principal Axes Alignment** - Molecules auto-rotated for optimal slicing orientation
- **xTB Integration** - Deformation density via xTB (Homebrew/conda), with install dialog
- **Detail Level Settings** - Simple mode (beginners) vs Advanced mode (shows π-bonds)
- **Fixed Contour Levels** - Advanced mode uses consistent levels across slices for comparability
- **i18n Ready** - Qt translation system in place (English only for now)

### What's Not Yet Complete

- Settings persistence
- App packaging (PyInstaller/cx_Freeze)

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
│   │   ├── cif_parser.py       # CIF file parsing
│   │   ├── cod_search.py       # COD + local molecule search
│   │   ├── density.py          # Electron density calculations
│   │   ├── contours.py         # Contour generation, PDF export
│   │   ├── generator.py        # Main generation pipeline
│   │   └── xtb_manager.py      # xTB download and execution
│   ├── ui/
│   │   ├── main_window.py      # Main application window
│   │   ├── molecule_viewer.py  # 2D structure viewer
│   │   └── slice_explorer.py   # 3D interactive slice preview
│   └── resources/
│       └── translations/       # i18n files
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
6. Unit tests

---

## Technical Notes

### Why xTB?

- Open source (LGPL-3.0), lightweight
- Handles organic molecules well
- Fast (seconds for small molecules)
- Auto-downloaded, not bundled (clean licensing)

### Deformation Density Calculation

The deformation density uses coordinates directly from the xTB cube file output to ensure perfect alignment between molecular and promolecule densities. The molecule is aligned to principal axes before xTB calculation so that z-slices cut parallel to the molecular plane (optimal for planar molecules like aromatics).

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

This tool gives students the endpoint of her process—density maps—and asks them to do what she did: interpret the contours as atoms and deduce chemical bonds.

---

## License

MIT (application code)

xTB is LGPL-3.0 (downloaded separately, not bundled)

---

*December 2024*

---

## Changelog

### v0.3.3 (December 2024)
- Fixed coordinate system alignment for custom plane selection
- Contours now properly align with 3D display when using 4-atom plane definition
- Unified right-handed coordinate convention (u × v = normal) across all modules
- Improved basis vector handling in density sampling and visualization

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
