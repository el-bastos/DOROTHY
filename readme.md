# CrystalSlice: Crystallography Teaching Tool

## Project Concept

Recreate the historical electron density contour maps from X-ray crystallography (the "museum plates" from the 1930s–1960s) as a hands-on teaching tool. Students physically stack transparent slices to visualize 3D electron density and deduce molecular structure—essentially reverse-engineering what pioneering crystallographers like Dorothy Hodgkin did.

### The Twist: Two Cubes

1. **Cube 1: Promolecule Density** — Superposition of spherical, non-interacting atoms at refined positions. Shows "here are the atoms."

2. **Cube 2: Deformation Density** — (Full molecular electron density) minus (promolecule density). Isolates bonding effects: buildup in bonding regions, lone pairs, depletion elsewhere.

This pairing teaches both structure determination and chemical bonding concepts.

---

## Technical Workflow

### Data Source

- **COD (Crystallography Open Database)** — Open access, no license required
- CIF files provide: unit cell parameters, space group, fractional coordinates, atom types
- User can search/browse COD and select molecules

### Calculation Pipeline

1. User selects molecule from COD
2. Backend extracts geometry from CIF
3. Orient molecule along principal axes
4. **xTB (GFN2-xTB)** single-point calculation → outputs electron density cube file
5. Calculate promolecule density (sum of tabulated spherical atomic densities)
6. Subtract promolecule from molecular density → deformation density
7. Slice both density cubes at user-defined z-levels (default: 15)
8. Generate PDF contour plots for printing

### Molecule Selection Criteria

- Orient molecule first (principal axes alignment)
- Only molecules that fit well in the slicing scheme
- One full molecule per output (not unit cell fragments)
- Small organic molecules only (initially)
- 10–20 atoms: complex enough to be a puzzle, simple enough to resolve

### Why xTB?

- Open source (LGPL-3.0), lightweight
- Handles organic molecules well
- Fast enough for on-demand calculation (seconds for small molecules like aspirin)
- Genuinely "black-boxable" for students

### Grid Resolution Presets

| Setting | Grid Spacing | Use Case |
|---------|-------------|----------|
| Coarse | ~0.2 Å | Quick preview, very clean lines |
| Medium | ~0.1 Å | Standard teaching use (default) |
| Fine | ~0.05 Å | Publication quality, subtle features |

---

## Physical Specifications

### Format

- **Default 15 slices** per cube (adjustable)
- **A5 size** (148 × 210 mm)
- Material: Adhesive transparent film applied to 2 mm transparent PVC sheets
- User handles physical assembly

### Contour Styling

**Promolecule density (Cube 1):**
- Solid black lines
- Evenly spaced contour levels
- Consistent scale across all slices

**Deformation density (Cube 2):**

| Feature | Meaning | Line Style |
|---------|---------|------------|
| **Positive** | Electron accumulation due to bonding | Solid lines |
| **Negative** | Electron depletion (moved to bonds) | Dashed/dotted lines |
| **Zero** | Omitted (noisy, distracting) | — |

**Color options:**
- Black and white (solid/dashed)
- Color mode: blue for positive, red for negative

### What Positive/Negative Deformation Shows

**Positive (excess density)** — electrons accumulated here *because* of bonding:
- Between bonded atoms (shared electrons in σ bonds)
- Above/below aromatic rings (π clouds)
- At lone pair positions (oxygen, nitrogen)

**Negative (depleted density)** — electrons moved away to participate in bonding:
- Close to nuclei along bond axes
- Opposite side from lone pairs

---

## Pilot Molecules

| Molecule | Teaching Value |
|----------|----------------|
| **Aspirin** | Aromatic ring, carbonyl, ester, carboxylic acid, ~21 atoms |
| Benzene | Simplest aromatic baseline |
| Naphthalene | Fused rings, shared bonding |
| Urea | Hydrogen bonding, very compact |
| Paracetamol | Amide vs ester comparison to aspirin |
| Caffeine | Student interest, multiple ring features |

---

## Correspondence with Experimental Data

### What's Experimental

The CIF coordinates from COD come from real diffraction experiments—actual X-rays hitting actual crystals. Atom positions are refined against measured intensities.

### What's Calculated

Both density cubes are calculated using quantum mechanics (xTB) to predict electron distribution given the experimental nuclear positions.

### Experimental Equivalent

High-resolution X-ray crystallography can measure deformation density directly using:
- Very high resolution data (<0.5 Å)
- Low temperature, excellent crystals
- Multipole refinement (Hansen-Coppens formalism)

Experimental deformation densities match DFT calculations remarkably well.

### Pedagogical Framing

- **Cube 1 (promolecule)**: "What crystallographers traditionally assume—spherical atoms"
- **Cube 2 (deformation)**: "What we now know is really there—and can measure with careful experiments or predict with quantum mechanics"

---

## Historical Context: Dorothy Hodgkin's Penicillin Protocol (1942–1945)

### Her Workflow

1. **Crystals** — Grew multiple crystal forms (K, Rb, Na salts of penicillin)
2. **Diffraction data** — X-ray photographs, intensities estimated visually
3. **Patterson maps** — Used |F|² only (no phases), compared K and Rb salts to locate heavy atoms
4. **Phase estimation** — Isomorphous replacement using heavy atom positions
5. **Fourier synthesis** — Initially by hand (Beevers-Lipson strips), later IBM punch-card tabulators
6. **Contour maps on Perspex** — Drew electron density contours on transparent sheets and stacked them
7. **Model building** — Wire and cork models from density interpretation

### Comparison

| Aspect | Hodgkin 1945 | CrystalSlice |
|--------|--------------|--------------|
| Atomic positions | Solved from experiment | Given from COD (experimental) |
| Electron density | Measured (via phases) | Calculated (xTB) |
| Deformation density | Not accessible | Calculated |
| Physical medium | Hand-drawn on Perspex | Printed film on PVC |
| Student task | Was the research itself | Reverse-engineer the structure |

The teaching tool gives students the endpoint of Hodgkin's process (density maps) and asks them to do what she did at the final stage—interpret the blobs as atoms and draw bonds.

---

## Application Architecture

### Stack

```
Python application
├── UI: PyQt6 or PySide6 (native desktop)
├── Chemistry: ASE (atomic simulation environment)
├── Density: NumPy + cube file handling
├── Contours: Matplotlib
├── PDF output: Matplotlib direct
├── i18n: Qt Linguist (built into Qt)
└── xTB: auto-downloaded on first run
```

### Why This Stack

| Aspect | Benefit |
|--------|---------|
| Python | Scientific ecosystem, maintainable |
| PyQt/PySide | True native look, no browser, cross-platform |
| Auto-download xTB | Zero user setup, clean licensing |
| Qt i18n | Mature system, supports multiple languages |
| Single package | One executable per platform |

### xTB Integration

**First run dialog:**
```
┌────────────────────────────────────────────────┐
│                                                │
│   Welcome to CrystalSlice                      │
│                                                │
│   This app needs a calculation engine (xTB)   │
│   to generate electron density maps.          │
│                                                │
│   Download now? (~20 MB, one time only)       │
│                                                │
│   xTB is free software from the Grimme group  │
│   at University of Bonn.                      │
│                                                │
│          [Download]    [Cancel]               │
│                                                │
└────────────────────────────────────────────────┘
```

- Downloads from official GitHub releases
- Stores in user's app data folder
- User always gets latest version
- No license entanglement (not bundling, just automating download)

### Licensing Note

xTB is **LGPL-3.0** licensed. By auto-downloading rather than bundling:
- Not distributing xTB directly
- User downloads from official source
- Authors get accurate download counts
- Clean separation of projects

Credits shown in app:
> "Electron density calculations powered by xTB, developed by the Grimme group (University of Bonn). https://github.com/grimme-lab/xtb"

### Distribution

```
CrystalSlice-Windows.zip     → extract, run .exe
CrystalSlice-macOS.dmg       → drag to Applications
CrystalSlice-Linux.AppImage  → chmod +x, run
```

Built with PyInstaller or cx_Freeze.

---

## User Interface

### Screen 1: Home / Search

```
┌─────────────────────────────────────────────────────┐
│  CrystalSlice                        [≡] [?] [⚙]   │
├─────────────────────────────────────────────────────┤
│                                                     │
│  ┌───────────────────────────────────┐              │
│  │ 🔍 Search molecule...             │  [Search]   │
│  └───────────────────────────────────┘              │
│                                                     │
│  Examples:  aspirin · benzene · caffeine · urea    │
│                                                     │
│  ─────────────────────────────────────────────────  │
│                                                     │
│  Recent:                                            │
│  ├─ Aspirin (COD 2300212)           [→]            │
│  └─ Naphthalene (COD 1503955)       [→]            │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### Screen 2: Search Results

```
┌─────────────────────────────────────────────────────┐
│  ← Results for "aspirin"                            │
├─────────────────────────────────────────────────────┤
│                                                     │
│  ┌─────────┐  Acetylsalicylic acid                 │
│  │  [3D]   │  C₉H₈O₄ · 21 atoms                    │
│  │  view   │  COD 2300212                          │
│  └─────────┘  ────────────────────────── [Select]  │
│                                                     │
│  ┌─────────┐  Aspirin polymorph II                 │
│  │  [3D]   │  C₉H₈O₄ · 21 atoms                    │
│  │  view   │  COD 7050228                          │
│  └─────────┘  ────────────────────────── [Select]  │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### Screen 3: Molecule Preview + Settings

```
┌─────────────────────────────────────────────────────┐
│  ← Aspirin                                          │
├─────────────────────────────────────────────────────┤
│                                                     │
│  ┌─────────────────────┐  Formula: C₉H₈O₄          │
│  │                     │  Atoms: 21                 │
│  │     [3D viewer]     │  Space group: P2₁/c       │
│  │     rotate/zoom     │                           │
│  │                     │  Source: COD 2300212      │
│  └─────────────────────┘                           │
│                                                     │
│  ─────────── Settings ───────────────────────────  │
│                                                     │
│  Resolution:   ○ Coarse  ● Medium  ○ Fine          │
│                                                     │
│  Output:       ☑ Promolecule density               │
│                ☑ Deformation density               │
│                                                     │
│  Color mode:   ○ Black/white  ● Color (red/blue)   │
│                                                     │
│  Slices:       15  [−] [+]                         │
│                                                     │
│               [ Generate PDFs ]                    │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### Screen 4: Processing

```
┌─────────────────────────────────────────────────────┐
│  Generating...                                      │
├─────────────────────────────────────────────────────┤
│                                                     │
│                                                     │
│         ████████████░░░░░░░░  60%                  │
│                                                     │
│         Running xTB calculation...                 │
│                                                     │
│                                                     │
│                    [Cancel]                        │
│                                                     │
└─────────────────────────────────────────────────────┘
```

Progress steps:
1. Fetching structure ✓
2. Running xTB calculation...
3. Calculating promolecule
4. Computing deformation density
5. Generating slices
6. Creating PDFs

### Screen 5: Results / Export

```
┌─────────────────────────────────────────────────────┐
│  ✓ Complete                                         │
├─────────────────────────────────────────────────────┤
│                                                     │
│  ┌──────────────────┐  ┌──────────────────┐        │
│  │ [slice preview]  │  │ [slice preview]  │        │
│  │                  │  │                  │        │
│  │  Promolecule     │  │  Deformation     │        │
│  │  Slice 8/15      │  │  Slice 8/15      │        │
│  └──────────────────┘  └──────────────────┘        │
│         [◀ ▶]                 [◀ ▶]                │
│                                                     │
│  ─────────────────────────────────────────────────  │
│                                                     │
│  [ Save all PDFs... ]   [ Save ZIP... ]            │
│                                                     │
│  [ ← New molecule ]                                │
│                                                     │
└─────────────────────────────────────────────────────┘
```

Preview allows scrolling through all slices before export.

### Settings Screen (⚙)

- Language selection (multi-language support)
- Default resolution preset
- Default slice count
- Default color mode
- Contour levels (advanced)
- xTB path (if user wants custom installation)
- About / credits

---

## Output Package

For each molecule, user downloads:

```
aspirin_crystalslice/
├── promolecule/
│   ├── slice_01.pdf
│   ├── slice_02.pdf
│   └── ... (15 files)
├── deformation/
│   ├── slice_01.pdf
│   ├── slice_02.pdf
│   └── ... (15 files)
└── info.txt (molecule details, settings used)
```

Or as a single ZIP file.

---

## Next Steps

1. Set up Python project structure
2. Implement COD search/fetch
3. xTB auto-download mechanism
4. Density calculation pipeline
5. Contour generation and PDF export
6. PyQt UI implementation
7. i18n setup
8. Packaging for all platforms
9. Testing with pilot molecules

---

*Discussion date: December 2024*
