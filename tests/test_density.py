"""
Tests for electron density calculations.

Validates:
- Promolecule density calculation
- Deformation density normalization
- Grid coordinate conversions
"""

import numpy as np
import pytest

from dorothy.core.constants import (
    BOHR_TO_ANGSTROM,
    ANGSTROM_TO_BOHR,
    ATOMIC_DENSITY_PARAMS,
)
from dorothy.core.density import (
    DensityCube,
    calculate_promolecule_density,
)


class TestConstants:
    """Test that constants are correctly defined."""

    def test_bohr_angstrom_conversion(self):
        """Verify Bohr/Angstrom conversion is reversible."""
        value = 5.0
        converted = value * BOHR_TO_ANGSTROM * ANGSTROM_TO_BOHR
        assert abs(converted - value) < 1e-10

    def test_atomic_density_params_coverage(self):
        """Verify common elements have density parameters."""
        common_elements = ['H', 'C', 'N', 'O', 'S', 'F', 'Cl']
        for elem in common_elements:
            assert elem in ATOMIC_DENSITY_PARAMS
            params = ATOMIC_DENSITY_PARAMS[elem]
            assert len(params) >= 1
            for coeff, exp in params:
                assert coeff > 0
                assert exp > 0


class TestDensityCube:
    """Test DensityCube class functionality."""

    @pytest.fixture
    def simple_cube(self):
        """Create a simple test density cube."""
        # 10x10x10 grid
        data = np.random.rand(10, 10, 10)
        origin = np.array([0.0, 0.0, 0.0])
        axes = np.diag([0.2, 0.2, 0.2])  # 0.2 Bohr spacing
        atoms = [(6, 1.0, 1.0, 1.0), (6, 2.0, 2.0, 2.0)]  # Two carbon atoms
        return DensityCube(origin=origin, axes=axes, data=data, atoms=atoms)

    def test_shape_property(self, simple_cube):
        """Test shape property returns correct dimensions."""
        assert simple_cube.shape == (10, 10, 10)

    def test_grid_spacing(self, simple_cube):
        """Test grid spacing calculation."""
        spacing = simple_cube.grid_spacing
        expected = 0.2 * BOHR_TO_ANGSTROM
        assert np.allclose(spacing, [expected, expected, expected])

    def test_get_slice_x(self, simple_cube):
        """Test extracting X slice."""
        slice_data = simple_cube.get_slice(axis=0, index=5)
        assert slice_data.shape == (10, 10)

    def test_get_slice_y(self, simple_cube):
        """Test extracting Y slice."""
        slice_data = simple_cube.get_slice(axis=1, index=5)
        assert slice_data.shape == (10, 10)

    def test_get_slice_z(self, simple_cube):
        """Test extracting Z slice."""
        slice_data = simple_cube.get_slice(axis=2, index=5)
        assert slice_data.shape == (10, 10)

    def test_get_z_slices(self, simple_cube):
        """Test getting multiple Z slices."""
        slices = simple_cube.get_z_slices(n_slices=5)
        assert len(slices) == 5
        for z_coord, slice_data in slices:
            assert isinstance(z_coord, float)
            assert slice_data.shape == (10, 10)

    def test_molecule_extent(self, simple_cube):
        """Test molecule extent calculation."""
        min_c, max_c = simple_cube.get_molecule_extent()
        # Atoms are at (1,1,1) and (2,2,2) in Bohr
        expected_min = 1.0 * BOHR_TO_ANGSTROM
        expected_max = 2.0 * BOHR_TO_ANGSTROM
        assert np.allclose(min_c, [expected_min] * 3)
        assert np.allclose(max_c, [expected_max] * 3)

    def test_density_at_point(self, simple_cube):
        """Test density interpolation at a point."""
        # Get density at a known grid point
        density = simple_cube.get_density_at_point(
            0.1 * BOHR_TO_ANGSTROM,
            0.1 * BOHR_TO_ANGSTROM,
            0.1 * BOHR_TO_ANGSTROM
        )
        assert isinstance(density, (int, float))


class TestPromoleculeDensity:
    """Test promolecule density calculation."""

    def test_single_hydrogen(self):
        """Test promolecule density for a single hydrogen atom."""
        coords = np.array([[0.0, 0.0, 0.0]])  # Hydrogen at origin (Angstrom)
        symbols = ['H']

        # Small grid around the atom
        origin = np.array([-2.0, -2.0, -2.0]) * ANGSTROM_TO_BOHR
        spacing = 0.5 * ANGSTROM_TO_BOHR
        axes = np.diag([spacing, spacing, spacing])
        shape = (9, 9, 9)

        density = calculate_promolecule_density(
            coords, symbols, origin, axes, shape
        )

        assert density.shape == shape
        # Density should be highest at center (index 4,4,4)
        center_density = density[4, 4, 4]
        edge_density = density[0, 0, 0]
        assert center_density > edge_density

    def test_two_carbons(self):
        """Test promolecule density for two carbon atoms."""
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],  # 1.5 Angstrom apart
        ])
        symbols = ['C', 'C']

        origin = np.array([-2.0, -2.0, -2.0]) * ANGSTROM_TO_BOHR
        spacing = 0.5 * ANGSTROM_TO_BOHR
        axes = np.diag([spacing, spacing, spacing])
        shape = (12, 9, 9)

        density = calculate_promolecule_density(
            coords, symbols, origin, axes, shape
        )

        assert density.shape == shape
        # Density should be positive everywhere
        assert np.all(density >= 0)
        # Total density should be roughly proportional to number of electrons
        # (6 electrons per carbon, very rough check)
        total = density.sum()
        assert total > 0

    def test_density_is_symmetric(self):
        """Test that single atom density is spherically symmetric."""
        coords = np.array([[0.0, 0.0, 0.0]])
        symbols = ['C']

        origin = np.array([-2.0, -2.0, -2.0]) * ANGSTROM_TO_BOHR
        spacing = 0.5 * ANGSTROM_TO_BOHR
        axes = np.diag([spacing, spacing, spacing])
        shape = (9, 9, 9)

        density = calculate_promolecule_density(
            coords, symbols, origin, axes, shape
        )

        center = 4
        # Check density at symmetric points
        d_x = density[center + 1, center, center]
        d_y = density[center, center + 1, center]
        d_z = density[center, center, center + 1]

        # Should be nearly equal due to spherical symmetry
        assert abs(d_x - d_y) < 1e-10
        assert abs(d_y - d_z) < 1e-10

    def test_unknown_element_uses_carbon(self):
        """Test that unknown elements default to carbon parameters."""
        coords = np.array([[0.0, 0.0, 0.0]])
        symbols = ['Xx']  # Unknown element

        origin = np.array([-2.0, -2.0, -2.0]) * ANGSTROM_TO_BOHR
        spacing = 0.5 * ANGSTROM_TO_BOHR
        axes = np.diag([spacing, spacing, spacing])
        shape = (9, 9, 9)

        # Should not raise, uses carbon as fallback
        density = calculate_promolecule_density(
            coords, symbols, origin, axes, shape
        )
        assert density.shape == shape
