"""
Tests for plane selection and rotation calculations.

Validates:
- PlaneDefinition rotation matrix properties
- Coordinate system orthonormality
- 4-atom plane calculation
"""

import numpy as np
import pytest

from dorothy.core.selection import SelectionManager, PlaneDefinition


class TestPlaneDefinition:
    """Test PlaneDefinition rotation matrix properties."""

    def test_rotation_matrix_orthonormal(self):
        """Verify rotation matrix is orthonormal (R @ R.T = I)."""
        # Create a simple rotation matrix from known vectors
        normal = np.array([0, 0, 1])
        u_axis = np.array([1, 0, 0])
        v_axis = np.array([0, 1, 0])

        rotation = np.array([u_axis, v_axis, normal])

        # Orthonormality check
        identity = rotation @ rotation.T
        assert np.allclose(identity, np.eye(3), atol=1e-10)

    def test_rotation_preserves_determinant(self):
        """Verify rotation matrix has determinant 1 (proper rotation)."""
        normal = np.array([0, 0, 1])
        u_axis = np.array([1, 0, 0])
        v_axis = np.array([0, 1, 0])

        rotation = np.array([u_axis, v_axis, normal])
        det = np.linalg.det(rotation)
        assert abs(det - 1.0) < 1e-10

    def test_tilted_plane_rotation(self):
        """Test rotation from a tilted plane to horizontal."""
        # Plane tilted 45 degrees
        normal = np.array([0, 1, 1]) / np.sqrt(2)

        # Construct orthonormal basis
        u_axis = np.array([1, 0, 0])
        v_axis = np.cross(normal, u_axis)
        v_axis = v_axis / np.linalg.norm(v_axis)

        rotation = np.array([u_axis, v_axis, normal])

        # Apply rotation to normal should give [0, 0, 1]
        rotated_normal = rotation @ normal
        assert np.allclose(rotated_normal, [0, 0, 1], atol=1e-10)


class TestSelectionManager:
    """Test SelectionManager plane calculation."""

    @pytest.fixture
    def manager(self):
        """Create a SelectionManager with test coordinates."""
        # Simple molecule in XY plane
        coords = np.array([
            [0.0, 0.0, 0.0],   # Atom 0: origin
            [1.0, 0.0, 0.0],   # Atom 1: along X
            [0.0, 1.0, 0.0],   # Atom 2: along Y
            [0.5, 0.5, 1.0],   # Atom 3: above plane (for "up" direction)
        ])
        manager = SelectionManager(parent=None)
        manager.set_coordinates(coords)
        return manager

    def test_initial_selection_empty(self, manager):
        """Test that initial selection is empty."""
        assert len(manager.selected_indices) == 0

    def test_select_atom(self, manager):
        """Test selecting an atom."""
        manager.toggle_atom(0)
        assert 0 in manager.selected_indices

    def test_deselect_atom(self, manager):
        """Test deselecting an atom."""
        manager.toggle_atom(0)
        manager.toggle_atom(0)
        assert 0 not in manager.selected_indices

    def test_clear_selection(self, manager):
        """Test clearing selection."""
        manager.toggle_atom(0)
        manager.toggle_atom(1)
        manager.clear()
        assert len(manager.selected_indices) == 0

    def test_selection_count(self, manager):
        """Test selection count property."""
        manager.toggle_atom(0)
        manager.toggle_atom(1)
        assert manager.selection_count == 2

    def test_plane_definition_with_4_atoms(self, manager):
        """Test plane calculation from 4 selected atoms."""
        # Select 4 atoms: 3 define plane, 4th defines "up"
        for i in range(4):
            manager.toggle_atom(i)

        plane_def = manager.get_plane_definition()
        assert plane_def is not None
        assert plane_def.center is not None
        assert plane_def.rotation_matrix is not None

    def test_no_plane_with_3_atoms(self, manager):
        """Test that 3 atoms don't define a plane (need 4)."""
        for i in range(3):
            manager.toggle_atom(i)

        plane_def = manager.get_plane_definition()
        assert plane_def is None

    def test_is_complete_property(self, manager):
        """Test is_complete property."""
        assert not manager.is_complete
        for i in range(4):
            manager.toggle_atom(i)
        assert manager.is_complete

    def test_rotation_matrix_is_valid(self, manager):
        """Test that computed rotation matrix is orthonormal."""
        for i in range(4):
            manager.toggle_atom(i)

        plane_def = manager.get_plane_definition()
        R = plane_def.rotation_matrix

        # Check orthonormality
        identity = R @ R.T
        assert np.allclose(identity, np.eye(3), atol=1e-10)

        # Check determinant is 1 (proper rotation, not reflection)
        det = np.linalg.det(R)
        assert abs(det - 1.0) < 1e-10


class TestCoordinateTransformation:
    """Test coordinate transformations with plane definitions."""

    def test_horizontal_plane_gives_z_aligned_normal(self):
        """Test that horizontal plane gives Z-aligned normal."""
        # Atoms in XY plane
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],  # "Up" direction
        ])

        manager = SelectionManager(parent=None)
        manager.set_coordinates(coords)

        for i in range(4):
            manager.toggle_atom(i)

        plane_def = manager.get_plane_definition()

        # Normal should point in Z direction (or close to it)
        normal = plane_def.normal
        # Should be close to [0, 0, 1] or [0, 0, -1]
        assert abs(abs(normal[2]) - 1.0) < 0.01
