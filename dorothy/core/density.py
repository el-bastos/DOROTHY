"""
Electron density calculations and cube file handling.

Handles:
- Parsing Gaussian cube files from xTB
- Calculating promolecule density (spherical atomic densities)
- Computing deformation density (molecular - promolecule)
"""

import numpy as np
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from dorothy.core.cif_parser import MoleculeStructure


# Atomic electron densities - simplified Slater-type functions
# Parameters: (coefficient, exponent) for sum of Gaussians approximation
# These give reasonable spherical atomic densities for visualization
ATOMIC_DENSITY_PARAMS = {
    'H': [(0.2829, 1.0), (0.4858, 0.3)],
    'C': [(2.0, 2.5), (4.0, 0.8)],
    'N': [(2.5, 2.8), (4.5, 0.9)],
    'O': [(3.0, 3.2), (5.0, 1.0)],
    'S': [(6.0, 2.0), (10.0, 0.5)],
    'P': [(5.0, 1.8), (10.0, 0.5)],
    'F': [(3.5, 3.5), (5.5, 1.1)],
    'Cl': [(7.0, 2.2), (10.0, 0.6)],
    'Br': [(15.0, 1.8), (20.0, 0.4)],
    'I': [(25.0, 1.5), (30.0, 0.35)],
}


@dataclass
class DensityCube:
    """3D electron density on a regular grid."""
    origin: np.ndarray  # (3,) origin in Bohr
    axes: np.ndarray  # (3, 3) axis vectors in Bohr
    data: np.ndarray  # (nx, ny, nz) density values
    atoms: list  # [(Z, x, y, z), ...] atomic numbers and positions

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.data.shape

    @property
    def grid_spacing(self) -> np.ndarray:
        """Grid spacing along each axis in Angstrom."""
        return np.linalg.norm(self.axes, axis=1) * 0.529177  # Bohr to Angstrom

    def get_slice(self, axis: int, index: int) -> np.ndarray:
        """Get a 2D slice perpendicular to given axis."""
        if axis == 0:
            return self.data[index, :, :]
        elif axis == 1:
            return self.data[:, index, :]
        else:
            return self.data[:, :, index]

    def get_z_slices(self, n_slices: int) -> list[tuple[float, np.ndarray]]:
        """
        Get evenly spaced slices along z-axis.

        Returns:
            List of (z_coordinate, 2D_array) tuples
        """
        nz = self.shape[2]
        indices = np.linspace(0, nz - 1, n_slices, dtype=int)

        slices = []
        for idx in indices:
            z_coord = self.origin[2] + idx * self.axes[2, 2]
            z_angstrom = z_coord * 0.529177  # Bohr to Angstrom
            slice_data = self.data[:, :, idx]
            slices.append((z_angstrom, slice_data))

        return slices


def parse_cube_file(filepath: Path) -> Optional[DensityCube]:
    """
    Parse a Gaussian cube file.

    Cube file format:
    - Line 1-2: Comments
    - Line 3: N_atoms, origin_x, origin_y, origin_z
    - Line 4-6: N_points, axis_vector for each axis
    - Next N_atoms lines: atomic_number, charge, x, y, z
    - Rest: density data
    """
    filepath = Path(filepath)
    if not filepath.exists():
        return None

    try:
        with open(filepath, 'r') as f:
            # Skip comments
            f.readline()
            f.readline()

            # Read header
            parts = f.readline().split()
            n_atoms = int(parts[0])
            origin = np.array([float(x) for x in parts[1:4]])

            axes = np.zeros((3, 3))
            n_points = []
            for i in range(3):
                parts = f.readline().split()
                n_points.append(int(parts[0]))
                axes[i] = [float(x) for x in parts[1:4]]

            # Read atoms
            atoms = []
            for _ in range(abs(n_atoms)):
                parts = f.readline().split()
                z = int(parts[0])
                pos = [float(x) for x in parts[2:5]]
                atoms.append((z, *pos))

            # Read density data
            data_flat = []
            for line in f:
                data_flat.extend([float(x) for x in line.split()])

            data = np.array(data_flat).reshape(n_points)

            return DensityCube(
                origin=origin,
                axes=axes,
                data=data,
                atoms=atoms,
            )

    except Exception as e:
        print(f"Error parsing cube file: {e}")
        return None


def calculate_promolecule_density(
    structure: MoleculeStructure,
    grid_origin: np.ndarray,
    grid_axes: np.ndarray,
    grid_shape: tuple[int, int, int],
) -> np.ndarray:
    """
    Calculate promolecule density (sum of spherical atomic densities).

    Args:
        structure: Molecule structure
        grid_origin: Origin of the grid (Bohr)
        grid_axes: Axis vectors (Bohr)
        grid_shape: Number of grid points (nx, ny, nz)

    Returns:
        3D array of promolecule density
    """
    coords = structure.get_cartesian_coords()  # Angstrom
    symbols = structure.get_symbols()

    # Convert to Bohr
    coords_bohr = coords / 0.529177

    nx, ny, nz = grid_shape
    density = np.zeros(grid_shape)

    # Create grid coordinates
    x = np.arange(nx) * grid_axes[0, 0] + grid_origin[0]
    y = np.arange(ny) * grid_axes[1, 1] + grid_origin[1]
    z = np.arange(nz) * grid_axes[2, 2] + grid_origin[2]

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Add contribution from each atom
    for sym, (ax, ay, az) in zip(symbols, coords_bohr):
        params = ATOMIC_DENSITY_PARAMS.get(sym, ATOMIC_DENSITY_PARAMS['C'])

        # Distance from this atom
        r_sq = (X - ax)**2 + (Y - ay)**2 + (Z - az)**2

        # Sum of Gaussians
        for coeff, exp in params:
            density += coeff * np.exp(-exp * r_sq)

    return density


def calculate_deformation_density(
    molecular_density: DensityCube,
    structure: MoleculeStructure,
) -> DensityCube:
    """
    Calculate deformation density = molecular - promolecule.

    Args:
        molecular_density: Density cube from xTB
        structure: Molecule structure

    Returns:
        Deformation density cube
    """
    promolecule = calculate_promolecule_density(
        structure,
        molecular_density.origin,
        molecular_density.axes,
        molecular_density.shape,
    )

    # Normalize promolecule to have same total electrons
    mol_sum = np.sum(molecular_density.data)
    pro_sum = np.sum(promolecule)
    if pro_sum > 0:
        promolecule *= mol_sum / pro_sum

    deformation_data = molecular_density.data - promolecule

    return DensityCube(
        origin=molecular_density.origin,
        axes=molecular_density.axes,
        data=deformation_data,
        atoms=molecular_density.atoms,
    )


def create_density_cube_from_structure(
    structure: MoleculeStructure,
    resolution: str = "medium",
) -> DensityCube:
    """
    Create a promolecule density cube directly from structure.

    Useful when xTB is not available - gives at least promolecule density.

    Args:
        structure: Molecule structure
        resolution: "coarse", "medium", or "fine"

    Returns:
        DensityCube with promolecule density
    """
    coords = structure.get_cartesian_coords()

    # Determine grid bounds (add padding)
    padding = 3.0  # Angstrom
    min_coords = coords.min(axis=0) - padding
    max_coords = coords.max(axis=0) + padding

    # Grid spacing based on resolution
    spacing_map = {"coarse": 0.2, "medium": 0.1, "fine": 0.05}
    spacing = spacing_map.get(resolution, 0.1)

    # Calculate grid dimensions
    extent = max_coords - min_coords
    n_points = (extent / spacing).astype(int) + 1

    # Convert to Bohr
    origin_bohr = min_coords / 0.529177
    spacing_bohr = spacing / 0.529177

    axes = np.diag([spacing_bohr, spacing_bohr, spacing_bohr])

    # Calculate promolecule density
    density = calculate_promolecule_density(
        structure,
        origin_bohr,
        axes,
        tuple(n_points),
    )

    # Create atoms list
    symbols = structure.get_symbols()
    coords_bohr = coords / 0.529177
    element_to_z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'S': 16, 'P': 15, 'F': 9, 'Cl': 17, 'Br': 35, 'I': 53}
    atoms = [(element_to_z.get(s, 6), *c) for s, c in zip(symbols, coords_bohr)]

    return DensityCube(
        origin=origin_bohr,
        axes=axes,
        data=density,
        atoms=atoms,
    )
