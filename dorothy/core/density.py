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

    def get_oriented_slices(
        self,
        n_slices: int,
        plane_definition
    ) -> tuple[list[tuple[float, np.ndarray]], dict]:
        """
        Get slices perpendicular to a custom plane.

        This method samples the density cube along slices that are
        perpendicular to the given plane normal, using the provided
        u/v axes for orientation and center for positioning.

        The coordinate system is:
        - normal: direction of slice progression (perpendicular to slice planes)
        - u: "right" direction in the slice plane (first axis of 2D grid)
        - v: "up" direction in the slice plane (second axis of 2D grid)

        These form a right-handed coordinate system: u × v = normal

        Args:
            n_slices: Number of slices to extract
            plane_definition: Either a PlaneDefinition object (from SelectionManager)
                            or a np.ndarray normal vector (backward compatible)

        Returns:
            Tuple of:
            - List of (offset, 2D_array) tuples where offset is the distance
              along the normal from the plane center
            - Dict with 'normal', 'u', 'v' basis vectors, 'center', 'half_extent', and 'n_points'
        """
        from scipy.interpolate import RegularGridInterpolator

        # Handle both PlaneDefinition objects and raw normal vectors
        if hasattr(plane_definition, 'normal'):
            # PlaneDefinition object - use provided basis vectors directly
            normal = np.asarray(plane_definition.normal, dtype=float)
            normal = normal / np.linalg.norm(normal)
            u = np.asarray(plane_definition.u_axis, dtype=float)
            u = u / np.linalg.norm(u)
            v = np.asarray(plane_definition.v_axis, dtype=float)
            v = v / np.linalg.norm(v)
            user_center = np.asarray(plane_definition.center, dtype=float)
            use_user_center = True
        else:
            # Raw normal vector (backward compatibility)
            normal = np.asarray(plane_definition, dtype=float)
            normal = normal / np.linalg.norm(normal)
            u = None
            v = None
            user_center = None
            use_user_center = False

        # Convert grid parameters to Angstrom
        bohr_to_ang = 0.529177
        origin_ang = self.origin * bohr_to_ang
        axes_ang = self.axes * bohr_to_ang

        nx, ny, nz = self.shape

        # Create coordinate arrays for the grid
        x_coords = origin_ang[0] + np.arange(nx) * axes_ang[0, 0]
        y_coords = origin_ang[1] + np.arange(ny) * axes_ang[1, 1]
        z_coords = origin_ang[2] + np.arange(nz) * axes_ang[2, 2]

        # Create interpolator
        interpolator = RegularGridInterpolator(
            (x_coords, y_coords, z_coords),
            self.data,
            method='linear',
            bounds_error=False,
            fill_value=0.0
        )

        # Determine center point
        if use_user_center:
            center_3d = user_center.copy()
        elif self.atoms:
            atom_coords = np.array([(a[1], a[2], a[3]) for a in self.atoms]) * bohr_to_ang
            center_3d = atom_coords.mean(axis=0)
        else:
            center_3d = np.array([
                (x_coords[0] + x_coords[-1]) / 2,
                (y_coords[0] + y_coords[-1]) / 2,
                (z_coords[0] + z_coords[-1]) / 2,
            ])

        # Calculate extent along normal for slice distribution
        if self.atoms:
            atom_coords = np.array([(a[1], a[2], a[3]) for a in self.atoms]) * bohr_to_ang
            projections = (atom_coords - center_3d) @ normal
            min_proj = projections.min() - 1.0  # Add padding
            max_proj = projections.max() + 1.0
        else:
            corners = np.array([
                [x_coords[0], y_coords[0], z_coords[0]],
                [x_coords[-1], y_coords[0], z_coords[0]],
                [x_coords[0], y_coords[-1], z_coords[0]],
                [x_coords[0], y_coords[0], z_coords[-1]],
                [x_coords[-1], y_coords[-1], z_coords[-1]],
            ])
            projections = (corners - center_3d) @ normal
            min_proj = projections.min()
            max_proj = projections.max()

        # Create orthonormal basis for the slice plane if not provided
        # Use consistent convention: u × v = normal (right-handed)
        if u is None or v is None:
            # Find a vector not parallel to normal
            if abs(normal[2]) < 0.9:
                ref = np.array([0.0, 0.0, 1.0])
            else:
                ref = np.array([1.0, 0.0, 0.0])

            # Create right-handed basis: u × v = normal
            v = np.cross(normal, ref)
            v = v / np.linalg.norm(v)
            u = np.cross(v, normal)
            u = u / np.linalg.norm(u)

        # Determine slice sampling resolution
        grid_spacing = min(
            abs(axes_ang[0, 0]),
            abs(axes_ang[1, 1]),
            abs(axes_ang[2, 2])
        )
        extent = max(
            x_coords[-1] - x_coords[0],
            y_coords[-1] - y_coords[0],
            z_coords[-1] - z_coords[0]
        )
        n_points = int(extent / grid_spacing) + 1
        n_points = min(n_points, 200)  # Cap for performance
        half_extent = extent / 2

        # Generate slice offsets (relative to center)
        offsets = np.linspace(min_proj, max_proj, n_slices)

        slices = []
        for offset in offsets:
            # Create sampling grid for this slice
            u_range = np.linspace(-half_extent, half_extent, n_points)
            v_range = np.linspace(-half_extent, half_extent, n_points)
            U, V = np.meshgrid(u_range, v_range, indexing='ij')

            # Calculate 3D coordinates for each sample point
            # Point at grid[i,j] = center + offset*normal + u_range[i]*u + v_range[j]*v
            slice_center = center_3d + normal * offset
            X = slice_center[0] + U * u[0] + V * v[0]
            Y = slice_center[1] + U * u[1] + V * v[1]
            Z = slice_center[2] + U * u[2] + V * v[2]

            # Sample the density
            points = np.stack([X, Y, Z], axis=-1)
            slice_data = interpolator(points)

            slices.append((offset, slice_data))

        # Return basis info for positioning slices in 3D
        # Include n_points to ensure slice_explorer uses the same grid
        basis_info = {
            'normal': normal,
            'u': u,
            'v': v,
            'center': center_3d,
            'extent': extent,
            'half_extent': half_extent,
            'n_points': n_points,
        }

        return slices, basis_info


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
) -> DensityCube:
    """
    Calculate deformation density = molecular - promolecule.

    Args:
        molecular_density: Density cube from xTB
        structure: Molecule structure (used only for element symbols)

    Returns:
        Deformation density cube
    """
    # Use atom coordinates from the cube file - these match exactly what xTB used
    # The cube file stores atoms as (Z, x, y, z) in Bohr
    z_to_symbol = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I'}
    symbols = [z_to_symbol.get(atom[0], 'C') for atom in molecular_density.atoms]
    # Convert from Bohr to Angstrom
    coords = np.array([[atom[1], atom[2], atom[3]] for atom in molecular_density.atoms]) * 0.529177

    promolecule = calculate_promolecule_density(
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
    pro_integral = np.sum(promolecule) * dV

    if pro_integral > 0:
        # Scale promolecule so integrated density matches molecular
        promolecule *= mol_integral / pro_integral

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
    align_to_principal_axes: bool = True,
) -> DensityCube:
    """
    Create a promolecule density cube directly from structure.

    Useful when xTB is not available - gives at least promolecule density.

    Args:
        structure: Molecule structure
        resolution: "coarse", "medium", or "fine"
        align_to_principal_axes: If True, rotate molecule so Z-slicing
            cuts parallel to the molecular plane (recommended for visualization)

    Returns:
        DensityCube with promolecule density
    """
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
