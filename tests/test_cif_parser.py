"""
Tests for CIF file parsing.

Validates:
- CIF value parsing (handling uncertainties)
- Space group extraction
- Atom coordinate parsing
"""

import numpy as np
import pytest
from pathlib import Path
import tempfile

from dorothy.core.cif_parser import (
    parse_cif_value,
    MoleculeStructure,
    CellParameters,
    Atom,
    parse_cif,
)


class TestParseCifValue:
    """Test CIF value parsing utility."""

    def test_simple_float(self):
        """Test parsing a simple float."""
        assert parse_cif_value("1.234") == 1.234

    def test_float_with_uncertainty(self):
        """Test parsing float with uncertainty in parentheses."""
        # CIF format: 1.234(5) means 1.234 ± 0.005
        assert parse_cif_value("1.234(5)") == 1.234

    def test_negative_float(self):
        """Test parsing negative float."""
        assert parse_cif_value("-0.5") == -0.5

    def test_scientific_notation(self):
        """Test parsing scientific notation."""
        result = parse_cif_value("1.5e-3")
        assert abs(result - 0.0015) < 1e-10

    def test_question_mark(self):
        """Test that '?' returns 0.0."""
        assert parse_cif_value("?") == 0.0

    def test_empty_string(self):
        """Test that empty string returns 0.0."""
        assert parse_cif_value("") == 0.0

    def test_invalid_value(self):
        """Test that invalid value returns 0.0."""
        assert parse_cif_value("invalid") == 0.0

    def test_integer(self):
        """Test parsing integer."""
        assert parse_cif_value("42") == 42.0


class TestMoleculeStructure:
    """Test MoleculeStructure class."""

    @pytest.fixture
    def simple_structure(self):
        """Create a simple test structure (water-like)."""
        atoms = [
            Atom(label='O1', symbol='O', x=0.0, y=0.0, z=0.0),
            Atom(label='H1', symbol='H', x=0.5, y=0.5, z=0.0),
            Atom(label='H2', symbol='H', x=-0.5, y=0.5, z=0.0),
        ]
        cell = CellParameters(a=10.0, b=10.0, c=10.0, alpha=90.0, beta=90.0, gamma=90.0)
        return MoleculeStructure(
            name="Water",
            formula="H2O",
            cell=cell,
            space_group="P1",
            atoms=atoms,
        )

    def test_get_symbols(self, simple_structure):
        """Test getting element symbols."""
        symbols = simple_structure.get_symbols()
        assert symbols == ['O', 'H', 'H']

    def test_atom_count(self, simple_structure):
        """Test atom count property."""
        assert simple_structure.atom_count == 3

    def test_get_cartesian_coords(self, simple_structure):
        """Test getting Cartesian coordinates."""
        coords = simple_structure.get_cartesian_coords()
        assert coords.shape == (3, 3)

    def test_coords_alignment(self, simple_structure):
        """Test principal axes alignment shifts to positive coordinates."""
        coords_aligned = simple_structure.get_cartesian_coords(align_to_principal_axes=True)
        assert coords_aligned.shape == (3, 3)
        # After alignment, coordinates are shifted to be non-negative
        # (minimum coords become origin for grid generation)
        min_coords = coords_aligned.min(axis=0)
        assert np.allclose(min_coords, [0, 0, 0], atol=1e-10)


class TestCellParameters:
    """Test CellParameters class."""

    def test_cubic_cell_matrix(self):
        """Test Cartesian matrix for cubic cell."""
        cell = CellParameters(a=10.0, b=10.0, c=10.0, alpha=90.0, beta=90.0, gamma=90.0)
        matrix = cell.to_cartesian_matrix()

        # For cubic cell, matrix should be diagonal
        assert matrix.shape == (3, 3)
        assert abs(matrix[0, 0] - 10.0) < 0.01
        assert abs(matrix[1, 1] - 10.0) < 0.01
        assert abs(matrix[2, 2] - 10.0) < 0.01


class TestCifParsing:
    """Test full CIF file parsing."""

    @pytest.fixture
    def minimal_cif_content(self):
        """Create minimal valid CIF content."""
        return """
data_test
_cell_length_a 10.000
_cell_length_b 10.000
_cell_length_c 10.000
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'P 1'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.0 0.0 0.0
C2 C 0.1 0.0 0.0
"""

    def test_parse_minimal_cif(self, minimal_cif_content):
        """Test parsing minimal CIF content."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
            f.write(minimal_cif_content)
            f.flush()

            structure = parse_cif(Path(f.name))

            assert structure is not None
            assert len(structure.atoms) == 2
            assert structure.cell.a == 10.0
            assert structure.space_group == 'P 1'

    def test_parse_with_uncertainties(self):
        """Test parsing CIF with uncertainty values."""
        content = """
data_test
_cell_length_a 10.000(5)
_cell_length_b 10.000(5)
_cell_length_c 10.000(5)
_cell_angle_alpha 90.0(1)
_cell_angle_beta 90.0(1)
_cell_angle_gamma 90.0(1)
_symmetry_space_group_name_H-M 'P 21/c'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.1234(5) 0.5678(5) 0.9012(5)
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
            f.write(content)
            f.flush()

            structure = parse_cif(Path(f.name))

            assert structure is not None
            assert abs(structure.atoms[0].x - 0.1234) < 0.0001
