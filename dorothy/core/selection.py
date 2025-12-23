"""
Atom selection model for interactive slice orientation.

Manages selection of up to 3 atoms to define a custom slicing plane,
with signals for synchronization between 2D and 3D views.
"""
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from PyQt6.QtCore import QObject, pyqtSignal


@dataclass
class AtomSelection:
    """Represents selected atoms for plane definition."""

    indices: list[int] = field(default_factory=list)
    max_atoms: int = 3

    @property
    def is_complete(self) -> bool:
        """Check if enough atoms are selected to define a plane."""
        return len(self.indices) >= self.max_atoms

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
            True if selection is now complete (3 atoms selected)
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


class SelectionManager(QObject):
    """
    Manages atom selection state with signals for view synchronization.

    Signals:
        selection_changed: Emitted when selection changes, with list of indices
        plane_defined: Emitted when 3 atoms selected, with plane normal vector
    """

    selection_changed = pyqtSignal(list)  # List of selected atom indices
    plane_defined = pyqtSignal(object)  # np.ndarray plane normal vector

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
        """Check if plane is fully defined."""
        return self._selection.is_complete

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
        """Calculate plane normal from 3 selected atoms and emit signal."""
        if self._coords is None or not self._selection.is_complete:
            return

        # Get coordinates of selected atoms
        p1 = self._coords[self._selection.indices[0]]
        p2 = self._coords[self._selection.indices[1]]
        p3 = self._coords[self._selection.indices[2]]

        # Calculate two vectors in the plane
        v1 = p2 - p1
        v2 = p3 - p1

        # Calculate normal via cross product
        normal = np.cross(v1, v2)

        # Normalize
        norm = np.linalg.norm(normal)
        if norm < 1e-10:
            # Degenerate case: atoms are collinear
            # Fall back to Z-axis
            normal = np.array([0.0, 0.0, 1.0])
        else:
            normal = normal / norm

        # Ensure consistent orientation (positive Z component)
        if normal[2] < 0:
            normal = -normal

        self.plane_defined.emit(normal)

    def get_plane_normal(self) -> Optional[np.ndarray]:
        """
        Calculate and return plane normal without emitting signal.

        Returns:
            Plane normal vector if 3 atoms selected, None otherwise
        """
        if self._coords is None or not self._selection.is_complete:
            return None

        p1 = self._coords[self._selection.indices[0]]
        p2 = self._coords[self._selection.indices[1]]
        p3 = self._coords[self._selection.indices[2]]

        v1 = p2 - p1
        v2 = p3 - p1
        normal = np.cross(v1, v2)

        norm = np.linalg.norm(normal)
        if norm < 1e-10:
            return np.array([0.0, 0.0, 1.0])

        normal = normal / norm
        if normal[2] < 0:
            normal = -normal

        return normal
