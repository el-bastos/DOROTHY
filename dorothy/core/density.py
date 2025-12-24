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
        Get evenly spaced slices along z-axis, distributed over molecule extent.

        Returns:
            List of (z_coordinate, 2D_array) tuples
        """
        bohr_to_ang = 0.529177
        nz = self.shape[2]

        # Determine Z range based on atom positions (molecule extent)
        if self.atoms:
            atom_z = np.array([a[3] for a in self.atoms]) * bohr_to_ang
            z_min = atom_z.min() - 0.5  # Add 0.5 Angstrom padding
            z_max = atom_z.max() + 0.5
        else:
            # Fallback to full cube extent
            z_min = self.origin[2] * bohr_to_ang
            z_max = (self.origin[2] + (nz - 1) * self.axes[2, 2]) * bohr_to_ang

        # Generate Z coordinates for slices
        z_coords_ang = np.linspace(z_min, z_max, n_slices)

        slices = []
        for z_ang in z_coords_ang:
            # Convert to grid index
            z_bohr = z_ang / bohr_to_ang
            idx = int(round((z_bohr - self.origin[2]) / self.axes[2, 2]))
            idx = max(0, min(idx, nz - 1))  # Clamp to valid range

            slice_data = self.data[:, :, idx]
            slices.append((z_ang, slice_data))

        return slices

    def get_molecule_extent(self) -> tuple[np.ndarray, np.ndarray]:
        """Get the min/max coordinates of atoms in Angstrom.

        Returns:
            (min_coords, max_coords) as (3,) arrays in Angstrom
        """
        bohr_to_ang = 0.529177
        if not self.atoms:
            # Fallback to grid extent
            min_c = self.origin * bohr_to_ang
            max_c = (self.origin + np.array(self.shape) * np.diag(self.axes)) * bohr_to_ang
            return min_c, max_c

        coords = np.array([[a[1], a[2], a[3]] for a in self.atoms]) * bohr_to_ang
        return coords.min(axis=0), coords.max(axis=0)

    def get_density_at_point(self, x: float, y: float, z: float) -> float:
        """Get interpolated density at a point in Angstrom coordinates.

        Args:
            x, y, z: Coordinates in Angstrom

        Returns:
            Density value (in e/Bohr³)
        """
        bohr_to_ang = 0.529177
        # Convert to Bohr
        x_bohr = x / bohr_to_ang
        y_bohr = y / bohr_to_ang
        z_bohr = z / bohr_to_ang

        # Convert to grid indices
        ix = (x_bohr - self.origin[0]) / self.axes[0, 0]
        iy = (y_bohr - self.origin[1]) / self.axes[1, 1]
        iz = (z_bohr - self.origin[2]) / self.axes[2, 2]

        # Clamp to valid range
        ix = max(0, min(ix, self.shape[0] - 1))
        iy = max(0, min(iy, self.shape[1] - 1))
        iz = max(0, min(iz, self.shape[2] - 1))

        # Simple nearest-neighbor for now
        return self.data[int(round(ix)), int(round(iy)), int(round(iz))]


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
    coords: np.ndarray,
    symbols: list[str],
    grid_origin: np.ndarray,
    grid_axes: np.ndarray,
    grid_shape: tuple[int, int, int],
) -> np.ndarray:
    """
    Calculate promolecule density (sum of spherical atomic densities).

    Args:
        coords: Atomic coordinates in Angstrom (N x 3)
        symbols: Element symbols for each atom
        grid_origin: Origin of the grid (Bohr)
        grid_axes: Axis vectors (Bohr)
        grid_shape: Number of grid points (nx, ny, nz)

    Returns:
        3D array of promolecule density
    """

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
) -> tuple[DensityCube, DensityCube]:
    """
    Calculate deformation density = molecular - promolecule.

    Args:
        molecular_density: Density cube from xTB
        structure: Molecule structure (used only for element symbols)

    Returns:
        Tuple of (promolecule_cube, deformation_cube) - both on the same grid
        as the molecular density for consistent animation blending.
    """
    # Use atom coordinates from the cube file - these match exactly what xTB used
    # The cube file stores atoms as (Z, x, y, z) in Bohr
    z_to_symbol = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I'}
    symbols = [z_to_symbol.get(atom[0], 'C') for atom in molecular_density.atoms]
    # Convert from Bohr to Angstrom
    coords = np.array([[atom[1], atom[2], atom[3]] for atom in molecular_density.atoms]) * 0.529177

    promolecule_data = calculate_promolecule_density(
        coords,
        symbols,
        molecular_density.origin,
        molecular_density.axes,
        molecular_density.shape,
    )

    # Calculate volume element for proper integration
    # Volume = det(axes) for each grid cell (in Bohr^3)
    dV = abs(np.linalg.det(molecular_density.axes))

    # Normalize promolecule to have same integrated electron count
    # ∫ρ dV = Σ ρ_i * dV, so we compare sums directly (dV cancels)
    mol_integral = np.sum(molecular_density.data) * dV
    pro_integral = np.sum(promolecule_data) * dV

    if pro_integral > 0:
        # Scale promolecule so integrated density matches molecular
        promolecule_data *= mol_integral / pro_integral

    deformation_data = molecular_density.data - promolecule_data

    promolecule_cube = DensityCube(
        origin=molecular_density.origin,
        axes=molecular_density.axes,
        data=promolecule_data,
        atoms=molecular_density.atoms,
    )

    deformation_cube = DensityCube(
        origin=molecular_density.origin,
        axes=molecular_density.axes,
        data=deformation_data,
        atoms=molecular_density.atoms,
    )

    return promolecule_cube, deformation_cube


def create_density_cube_from_structure(
    structure: MoleculeStructure,
    resolution: str = "medium",
    align_to_principal_axes: bool = True,
    plane_definition=None,
) -> DensityCube:
    """
    Create a promolecule density cube directly from structure.

    Useful when xTB is not available - gives at least promolecule density.

    Args:
        structure: Molecule structure
        resolution: "coarse", "medium", or "fine"
        align_to_principal_axes: If True and no plane_definition, rotate molecule
            so Z-slicing cuts parallel to the molecular plane
        plane_definition: If provided, rotate coordinates so the user-selected plane
                         becomes horizontal (XY plane). Takes precedence over
                         align_to_principal_axes.

    Returns:
        DensityCube with promolecule density
    """
    if plane_definition is not None:
        # Use custom plane rotation: rotate so plane normal -> Z axis
        # IMPORTANT: The plane_definition was calculated from principal-axes-aligned
        # coordinates, so we must start from those same aligned coordinates
        coords = structure.get_cartesian_coords(align_to_principal_axes=True)
        center = plane_definition.center
        rotation = plane_definition.rotation_matrix
        # Center on plane center, then rotate
        coords_centered = coords - center
        coords = coords_centered @ rotation.T
    else:
        coords = structure.get_cartesian_coords(align_to_principal_axes=align_to_principal_axes)

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
    symbols = structure.get_symbols()
    density = calculate_promolecule_density(
        coords,
        symbols,
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
