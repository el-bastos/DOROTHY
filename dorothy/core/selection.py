"""
Atom selection model for interactive slice orientation.

Manages selection of up to 4 atoms to define a custom slicing plane:
- 3 atoms define the plane orientation
- 4th atom defines the "up" direction in the slice coordinate system
- Slice center is the centroid of all 4 atoms

With signals for synchronization between 2D and 3D views.
"""
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from PyQt6.QtCore import QObject, pyqtSignal


@dataclass
class AtomSelection:
    """Represents selected atoms for plane definition."""

    indices: list[int] = field(default_factory=list)
    max_atoms: int = 4  # 3 for plane + 1 for orientation

    @property
    def is_complete(self) -> bool:
        """Check if enough atoms are selected to define a plane with orientation."""
        return len(self.indices) >= self.max_atoms

    @property
    def can_define_plane(self) -> bool:
        """Check if at least 3 atoms are selected (minimum for a plane)."""
        return len(self.indices) >= 3

    @property
    def count(self) -> int:
        """Number of currently selected atoms."""
        return len(self.indices)

    def toggle(self, index: int) -> bool:
        """
        Toggle atom selection.

        Args:
            index: Atom index to toggle

        Returns:
            True if selection is now complete (4 atoms selected)
        """
        if index in self.indices:
            self.indices.remove(index)
            return False

        if len(self.indices) < self.max_atoms:
            self.indices.append(index)

        return self.is_complete

    def clear(self):
        """Clear all selections."""
        self.indices.clear()


@dataclass
class PlaneDefinition:
    """Complete plane definition from 4 atoms.

    The rotation_matrix transforms coordinates so that:
    - The selected plane becomes the XY plane
    - Z-axis aligns with the plane normal
    - Slicing along Z will cut parallel to the user-selected plane
    """
    normal: np.ndarray  # Unit normal to the plane
    u_axis: np.ndarray  # "Right" direction in the slice
    v_axis: np.ndarray  # "Up" direction in the slice
    center: np.ndarray  # Centroid of the 4 selected atoms
    rotation_matrix: np.ndarray = field(default_factory=lambda: np.eye(3))  # 3x3 rotation matrix


class SelectionManager(QObject):
    """
    Manages atom selection state with signals for view synchronization.

    4-atom selection:
    - Atoms 1-3 define the plane (via cross product)
    - Atom 4 defines the "up" direction in the slice
    - Center is the centroid of all 4 atoms

    Signals:
        selection_changed: Emitted when selection changes, with list of indices
        plane_defined: Emitted when 4 atoms selected, with PlaneDefinition
    """

    selection_changed = pyqtSignal(list)  # List of selected atom indices
    plane_defined = pyqtSignal(object)  # PlaneDefinition object

    def __init__(self, parent=None):
        super().__init__(parent)
        self._selection = AtomSelection()
        self._coords: Optional[np.ndarray] = None

    @property
    def selected_indices(self) -> list[int]:
        """Get current selection as list of indices."""
        return self._selection.indices.copy()

    @property
    def is_complete(self) -> bool:
        """Check if plane is fully defined (4 atoms)."""
        return self._selection.is_complete

    @property
    def selection_count(self) -> int:
        """Number of currently selected atoms."""
        return self._selection.count

    def set_coordinates(self, coords: np.ndarray):
        """
        Set the coordinate array for plane calculation.

        Args:
            coords: (N, 3) array of atom coordinates in Angstrom
        """
        self._coords = coords
        self.clear()

    def toggle_atom(self, index: int):
        """
        Toggle atom selection and emit appropriate signals.

        Args:
            index: Atom index to toggle
        """
        if self._coords is None:
            return

        if index < 0 or index >= len(self._coords):
            return

        complete = self._selection.toggle(index)
        self.selection_changed.emit(self._selection.indices.copy())

        if complete:
            self._calculate_and_emit_plane()

    def clear(self):
        """Clear selection and emit signal."""
        self._selection.clear()
        self.selection_changed.emit([])

    def _calculate_and_emit_plane(self):
        """Calculate plane definition from 4 selected atoms and emit signal."""
        if self._coords is None or not self._selection.is_complete:
            return

        plane_def = self._calculate_plane_definition()
        if plane_def:
            self.plane_defined.emit(plane_def)

    def _calculate_plane_definition(self) -> Optional[PlaneDefinition]:
        """
        Calculate complete plane definition from 4 atoms.

        - Atoms 1-3: Define the plane via cross product
        - Atom 4: Defines the "up" direction (v_axis points toward atom 4)
        - Center: Centroid of all 4 atoms

        The resulting coordinate system is right-handed: u × v = normal
        This ensures consistency with how density.py samples and
        slice_explorer.py displays the slices.
        """
        if self._coords is None or not self._selection.is_complete:
            return None

        # Get coordinates of selected atoms
        p1 = self._coords[self._selection.indices[0]]
        p2 = self._coords[self._selection.indices[1]]
        p3 = self._coords[self._selection.indices[2]]
        p4 = self._coords[self._selection.indices[3]]

        # Calculate plane normal from first 3 atoms
        v1 = p2 - p1
        v2 = p3 - p1
        normal = np.cross(v1, v2)

        norm = np.linalg.norm(normal)
        if norm < 1e-10:
            # Degenerate case: atoms are collinear
            normal = np.array([0.0, 0.0, 1.0])
        else:
            normal = normal / norm

        # Center is centroid of all 4 atoms
        center = (p1 + p2 + p3 + p4) / 4.0

        # Atom 4 defines the "up" direction (v_axis)
        # Project p4 onto the plane and use direction from center to projection
        p4_rel = p4 - center
        p4_in_plane = p4_rel - np.dot(p4_rel, normal) * normal

        if np.linalg.norm(p4_in_plane) < 1e-10:
            # p4 is directly above/below center, use arbitrary direction
            if abs(normal[2]) < 0.9:
                v_axis = np.array([0.0, 0.0, 1.0])
            else:
                v_axis = np.array([1.0, 0.0, 0.0])
            v_axis = v_axis - np.dot(v_axis, normal) * normal
            v_axis = v_axis / np.linalg.norm(v_axis)
        else:
            v_axis = p4_in_plane / np.linalg.norm(p4_in_plane)

        # u_axis completes the right-handed system: u × v = normal
        # This means u = v × normal (not normal × v)
        u_axis = np.cross(v_axis, normal)
        u_axis = u_axis / np.linalg.norm(u_axis)

        # Build rotation matrix that transforms coordinates so:
        # - u_axis -> X axis (1, 0, 0)
        # - v_axis -> Y axis (0, 1, 0)
        # - normal -> Z axis (0, 0, 1)
        # The rotation matrix R has rows = target basis vectors
        # coords_rotated = (coords - center) @ R.T
        rotation_matrix = np.array([u_axis, v_axis, normal])

        return PlaneDefinition(
            normal=normal,
            u_axis=u_axis,
            v_axis=v_axis,
            center=center,
            rotation_matrix=rotation_matrix,
        )

    def get_plane_definition(self) -> Optional[PlaneDefinition]:
        """
        Calculate and return plane definition without emitting signal.

        Returns:
            PlaneDefinition if 4 atoms selected, None otherwise
        """
        return self._calculate_plane_definition()

    def get_plane_normal(self) -> Optional[np.ndarray]:
        """
        Get just the plane normal (for backward compatibility).

        Returns:
            Plane normal vector if 4 atoms selected, None otherwise
        """
        plane_def = self._calculate_plane_definition()
        if plane_def:
            return plane_def.normal
        return None
