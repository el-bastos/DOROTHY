"""
Physical and application constants for Dorothy.

Centralizes magic numbers used throughout the codebase.
"""

# =============================================================================
# Physical Constants
# =============================================================================

# Bohr radius in Angstroms (a₀ = 0.529177 Å)
# Used for unit conversion between Bohr (atomic units) and Angstrom
BOHR_TO_ANGSTROM = 0.529177
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM

# =============================================================================
# Element Data
# =============================================================================

# Atomic number to symbol mapping
Z_TO_SYMBOL = {
    1: 'H', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
    14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
    33: 'As', 34: 'Se', 35: 'Br', 53: 'I',
}

# Symbol to atomic number mapping
SYMBOL_TO_Z = {v: k for k, v in Z_TO_SYMBOL.items()}

# Covalent radii in Angstroms (for bond detection)
# Source: Cordero et al., Dalton Trans., 2008, 2832-2838
COVALENT_RADII = {
    1: 0.31, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57,
    14: 1.11, 15: 1.07, 16: 1.05, 17: 1.02,
    33: 1.19, 34: 1.20, 35: 1.20, 53: 1.39,
}

# CPK colors for element visualization (hex)
ELEMENT_COLORS = {
    1: '#FFFFFF',   # H - white
    5: '#FFB5B5',   # B - salmon
    6: '#909090',   # C - gray
    7: '#3050F8',   # N - blue
    8: '#FF0D0D',   # O - red
    9: '#90E050',   # F - green
    14: '#F0C8A0',  # Si - tan
    15: '#FF8000',  # P - orange
    16: '#FFFF30',  # S - yellow
    17: '#1FF01F',  # Cl - green
    33: '#BD80E3',  # As - violet
    34: '#FFA100',  # Se - orange
    35: '#A62929',  # Br - brown
    53: '#940094',  # I - purple
}

# =============================================================================
# Density Calculation Parameters
# =============================================================================

# Grid spacing in Angstroms for different resolution levels
GRID_SPACING = {
    'coarse': 0.2,   # Quick preview
    'medium': 0.1,   # Default teaching use
    'fine': 0.05,    # Publication quality
}

# Atomic electron densities - simplified Slater-type functions
# Parameters: list of (coefficient, exponent) for sum of Gaussians approximation
# These give reasonable spherical atomic densities for visualization
# Note: These are approximate values tuned for visual representation,
# not rigorous quantum mechanical calculations
ATOMIC_DENSITY_PARAMS = {
    'H': [(0.2829, 1.0), (0.4858, 0.3)],
    'B': [(1.5, 2.0), (3.0, 0.7)],
    'C': [(2.0, 2.5), (4.0, 0.8)],
    'N': [(2.5, 2.8), (4.5, 0.9)],
    'O': [(3.0, 3.2), (5.0, 1.0)],
    'Si': [(5.0, 1.6), (9.0, 0.5)],
    'P': [(5.0, 1.8), (10.0, 0.5)],
    'S': [(6.0, 2.0), (10.0, 0.5)],
    'F': [(3.5, 3.5), (5.5, 1.1)],
    'Cl': [(7.0, 2.2), (10.0, 0.6)],
    'As': [(13.0, 1.8), (20.0, 0.5)],
    'Se': [(14.0, 1.8), (20.0, 0.5)],
    'Br': [(15.0, 1.8), (20.0, 0.4)],
    'I': [(25.0, 1.5), (30.0, 0.35)],
}

# =============================================================================
# Contour Level Settings
# =============================================================================

# Fixed contour levels for advanced mode (e/Bohr³)
# These levels are designed to show both σ and π bonding features
CONTOUR_LEVELS_PROMOLECULE = [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25]

# Deformation density contour levels (positive values; negative are mirrored)
# Lower thresholds capture diffuse π-density above/below molecular plane
CONTOUR_LEVELS_DEFORMATION = [0.01, 0.02, 0.04, 0.08, 0.12]

# Base levels for dynamic contour scaling
CONTOUR_LEVELS_BASE = [0.005, 0.01, 0.02, 0.04, 0.08]

# =============================================================================
# Visualization Settings
# =============================================================================

# Default number of slices for density visualization
DEFAULT_N_SLICES = 15

# Padding around molecule for grid bounds (Angstroms)
GRID_PADDING = 3.0

# Slice thickness for atom display (Angstroms)
# Atoms within this distance of a slice plane are shown
SLICE_ATOM_THICKNESS = 0.5

# Molecule extent padding (Angstroms)
# Added to min/max atom coordinates when determining slice range
MOLECULE_EXTENT_PADDING = 0.5

# =============================================================================
# Output Settings
# =============================================================================

# A5 page size in inches (148 × 210 mm)
PAGE_SIZE_A5 = (5.83, 8.27)
