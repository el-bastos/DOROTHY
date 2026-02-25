"""
2D Molecule viewer widget using matplotlib.

Supports interactive atom picking for plane definition and 3D-like rotation.
"""

import numpy as np
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QSizePolicy
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from dorothy.core.cif_parser import MoleculeStructure


# Element colors (CPK coloring scheme)
ELEMENT_COLORS = {
    'H': '#FFFFFF',
    'C': '#909090',
    'N': '#3050F8',
    'O': '#FF0D0D',
    'S': '#FFFF30',
    'P': '#FF8000',
    'F': '#90E050',
    'Cl': '#1FF01F',
    'Br': '#A62929',
    'I': '#940094',
}

# Covalent radii (Å) for bond detection
COVALENT_RADII = {
    'H': 0.31,
    'C': 0.76,
    'N': 0.71,
    'O': 0.66,
    'S': 1.05,
    'P': 1.07,
    'F': 0.57,
    'Cl': 1.02,
    'Br': 1.20,
    'I': 1.39,
}


class MoleculeCanvas(FigureCanvas):
    """Canvas for displaying 2D molecule structure with atom picking and rotation."""

    # Signal emitted when an atom is clicked (passes atom index)
    atom_picked = pyqtSignal(int)

    def __init__(self, parent=None):
        self.fig = Figure(figsize=(5, 5), dpi=100)
        self.fig.patch.set_facecolor('#f8f8f8')
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._structure = None
        self._selected_indices: list[int] = []
        self._pick_enabled = False
        self._atom_positions: list[tuple[float, float]] = []  # Store for hit testing

        # Rotation state (elevation and azimuth angles in degrees)
        self._elev = 0.0  # Rotation around X axis (tilt up/down)
        self._azim = 0.0  # Rotation around Z axis (spin left/right)

        # Mouse drag state
        self._dragging = False
        self._last_mouse_pos = None

        # Connect mouse events
        self.mpl_connect('button_press_event', self._on_click)
        self.mpl_connect('button_release_event', self._on_release)
        self.mpl_connect('motion_notify_event', self._on_motion)

        self._clear_plot()

    def _clear_plot(self):
        """Clear the plot."""
        self.ax.clear()
        self.ax.set_aspect('equal')
        self.ax.axis('off')
        self.fig.tight_layout()

    def _rotation_matrix(self, elev: float, azim: float) -> np.ndarray:
        """
        Create a rotation matrix for the given elevation and azimuth angles.

        Args:
            elev: Elevation angle in degrees (rotation around X axis)
            azim: Azimuth angle in degrees (rotation around Z axis)

        Returns:
            3x3 rotation matrix
        """
        elev_rad = np.radians(elev)
        azim_rad = np.radians(azim)

        # Rotation around Z axis (azimuth)
        Rz = np.array([
            [np.cos(azim_rad), -np.sin(azim_rad), 0],
            [np.sin(azim_rad), np.cos(azim_rad), 0],
            [0, 0, 1]
        ])

        # Rotation around X axis (elevation)
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(elev_rad), -np.sin(elev_rad)],
            [0, np.sin(elev_rad), np.cos(elev_rad)]
        ])

        # Combined rotation: first azimuth, then elevation
        return Rx @ Rz

    def set_structure(self, structure: MoleculeStructure):
        """Display a molecule structure."""
        self._structure = structure
        self._draw_molecule()

    def _draw_molecule(self):
        """Draw the molecule structure."""
        self._clear_plot()

        if not self._structure or not self._structure.atoms:
            self.ax.text(0.5, 0.5, 'No structure data',
                        ha='center', va='center', transform=self.ax.transAxes,
                        fontsize=12, color='#888')
            self.draw()
            return

        # Get Cartesian coordinates
        coords = self._structure.get_cartesian_coords()
        symbols = self._structure.get_symbols()

        if len(coords) == 0:
            return

        # Center the molecule first
        coords = coords - coords.mean(axis=0)

        # Apply rotation
        R = self._rotation_matrix(self._elev, self._azim)
        coords = coords @ R.T

        # Project to 2D (use x, y coordinates after rotation)
        x = coords[:, 0]
        y = coords[:, 1]

        # Find bonds based on distance (use original unrotated coords for distance calc)
        orig_coords = self._structure.get_cartesian_coords()
        bonds = self._find_bonds(orig_coords, symbols)

        # Sort atoms by Z coordinate (depth) for proper drawing order
        z = coords[:, 2]
        depth_order = np.argsort(z)  # Back to front

        # Draw bonds first (so they're behind atoms)
        for i, j in bonds:
            # Depth-based alpha for 3D effect
            avg_z = (z[i] + z[j]) / 2
            z_range = z.max() - z.min() if z.max() != z.min() else 1.0
            depth_factor = (avg_z - z.min()) / z_range if z_range > 0 else 0.5
            alpha = 0.4 + 0.6 * depth_factor  # Range 0.4 to 1.0
            linewidth = 1.5 + 1.0 * depth_factor  # Range 1.5 to 2.5
            self.ax.plot([x[i], x[j]], [y[i], y[j]],
                        color='#404040', linewidth=linewidth, alpha=alpha, zorder=1)

        # Store atom positions for hit testing
        self._atom_positions = list(zip(x, y))

        # Calculate depth factors for all atoms
        z_range = z.max() - z.min() if z.max() != z.min() else 1.0
        depth_factors = (z - z.min()) / z_range if z_range > 0 else np.full_like(z, 0.5)

        # Draw atoms in depth order (back to front)
        for idx in depth_order:
            xi, yi = x[idx], y[idx]
            sym = symbols[idx]
            depth_factor = depth_factors[idx]

            color = ELEMENT_COLORS.get(sym, '#808080')
            # Size based on element and depth (closer = larger)
            base_size = 150 if sym == 'H' else 300
            size = base_size * (0.6 + 0.8 * depth_factor)  # Range 0.6x to 1.4x

            # Highlight selected atoms
            if idx in self._selected_indices:
                edgecolor = '#FF0000'
                linewidth = 3
            else:
                edgecolor = '#404040'
                linewidth = 1

            # Alpha based on depth
            alpha = 0.5 + 0.5 * depth_factor  # Range 0.5 to 1.0

            # Draw atom circle
            self.ax.scatter(xi, yi, s=size, c=color, edgecolors=edgecolor,
                           linewidths=linewidth, zorder=2 + depth_factor, alpha=alpha)

            # Label (skip H for cleaner view)
            if sym != 'H':
                self.ax.annotate(sym, (xi, yi), ha='center', va='center',
                               fontsize=8, fontweight='bold', zorder=3 + depth_factor,
                               alpha=alpha)

        # Set axis limits with padding
        padding = 1.0
        self.ax.set_xlim(x.min() - padding, x.max() + padding)
        self.ax.set_ylim(y.min() - padding, y.max() + padding)

        self.fig.tight_layout()
        self.draw()

    def _find_bonds(self, coords: np.ndarray, symbols: list[str]) -> list[tuple[int, int]]:
        """Find bonds based on interatomic distances."""
        bonds = []
        n_atoms = len(symbols)

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                # Calculate distance
                dist = np.linalg.norm(coords[i] - coords[j])

                # Get covalent radii
                r1 = COVALENT_RADII.get(symbols[i], 1.5)
                r2 = COVALENT_RADII.get(symbols[j], 1.5)

                # Bond if distance < sum of covalent radii + tolerance
                max_bond_dist = (r1 + r2) * 1.3  # 30% tolerance
                if dist < max_bond_dist:
                    bonds.append((i, j))

        return bonds

    def _on_click(self, event):
        """Handle mouse click: left = atom picking, right = rotation drag."""
        if event.inaxes != self.ax:
            return

        # Left button: atom picking only
        if event.button == 1:
            if self._pick_enabled and self._atom_positions:
                click_x, click_y = event.x, event.y
                min_dist = float('inf')
                closest_idx = -1

                for idx, (ax, ay) in enumerate(self._atom_positions):
                    x_disp, y_disp = self.ax.transData.transform((ax, ay))
                    dist = np.sqrt((x_disp - click_x) ** 2 + (y_disp - click_y) ** 2)
                    if dist < min_dist:
                        min_dist = dist
                        closest_idx = idx

                if closest_idx >= 0 and min_dist < 25:
                    self.atom_picked.emit(closest_idx)
            return

        # Right button: start rotation drag
        if event.button == 3:
            self._dragging = True
            self._last_mouse_pos = (event.x, event.y)

    def _on_release(self, event):
        """Handle mouse button release."""
        self._dragging = False
        self._last_mouse_pos = None

    def _on_motion(self, event):
        """Handle mouse motion for rotation."""
        if not self._dragging or self._last_mouse_pos is None:
            return

        if event.x is None or event.y is None:
            return

        # Calculate mouse movement
        dx = event.x - self._last_mouse_pos[0]
        dy = event.y - self._last_mouse_pos[1]

        # Update angles (sensitivity factor)
        sensitivity = 0.5
        self._azim += dx * sensitivity
        self._elev -= dy * sensitivity  # Inverted for natural feel

        # Clamp elevation to avoid gimbal lock issues
        self._elev = max(-90, min(90, self._elev))

        # Normalize azimuth
        self._azim = self._azim % 360

        self._last_mouse_pos = (event.x, event.y)

        # Redraw
        if self._structure:
            self._draw_molecule()

    def set_selection(self, indices: list[int]):
        """
        Update visual selection.

        Args:
            indices: List of atom indices to highlight
        """
        self._selected_indices = indices.copy()
        if self._structure:
            self._draw_molecule()

    def set_pick_enabled(self, enabled: bool):
        """Enable or disable atom picking."""
        self._pick_enabled = enabled

    def get_selection(self) -> list[int]:
        """Get current selection."""
        return self._selected_indices.copy()

    def rotate_left(self):
        """Rotate view left by 15 degrees."""
        self._azim -= 15
        self._azim = self._azim % 360
        if self._structure:
            self._draw_molecule()

    def rotate_right(self):
        """Rotate view right by 15 degrees."""
        self._azim += 15
        self._azim = self._azim % 360
        if self._structure:
            self._draw_molecule()

    def rotate_up(self):
        """Rotate view up by 15 degrees."""
        self._elev = min(90, self._elev + 15)
        if self._structure:
            self._draw_molecule()

    def rotate_down(self):
        """Rotate view down by 15 degrees."""
        self._elev = max(-90, self._elev - 15)
        if self._structure:
            self._draw_molecule()

    def reset_rotation(self):
        """Reset rotation to default view."""
        self._elev = 0.0
        self._azim = 0.0
        if self._structure:
            self._draw_molecule()


class MoleculeViewer(QWidget):
    """Widget for viewing molecule structures with rotation controls."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(5)

        # Canvas
        self.canvas = MoleculeCanvas(self)
        layout.addWidget(self.canvas, stretch=1)

        # Rotation controls
        controls = QHBoxLayout()
        controls.setSpacing(5)

        # Add stretch to center the controls
        controls.addStretch()

        rotate_label = QPushButton("Rotate:")
        rotate_label.setFlat(True)
        rotate_label.setEnabled(False)
        rotate_label.setStyleSheet("font-weight: bold; font-size: 9pt; color: #666;")
        controls.addWidget(rotate_label)

        rotate_left_btn = QPushButton("\u2190")  # Left arrow
        rotate_left_btn.setMaximumWidth(30)
        rotate_left_btn.setToolTip("Rotate left (or drag mouse)")
        rotate_left_btn.clicked.connect(self.canvas.rotate_left)
        controls.addWidget(rotate_left_btn)

        rotate_right_btn = QPushButton("\u2192")  # Right arrow
        rotate_right_btn.setMaximumWidth(30)
        rotate_right_btn.setToolTip("Rotate right (or drag mouse)")
        rotate_right_btn.clicked.connect(self.canvas.rotate_right)
        controls.addWidget(rotate_right_btn)

        rotate_up_btn = QPushButton("\u2191")  # Up arrow
        rotate_up_btn.setMaximumWidth(30)
        rotate_up_btn.setToolTip("Rotate up (or drag mouse)")
        rotate_up_btn.clicked.connect(self.canvas.rotate_up)
        controls.addWidget(rotate_up_btn)

        rotate_down_btn = QPushButton("\u2193")  # Down arrow
        rotate_down_btn.setMaximumWidth(30)
        rotate_down_btn.setToolTip("Rotate down (or drag mouse)")
        rotate_down_btn.clicked.connect(self.canvas.rotate_down)
        controls.addWidget(rotate_down_btn)

        reset_btn = QPushButton("Reset")
        reset_btn.setMaximumWidth(50)
        reset_btn.setToolTip("Reset to default view")
        reset_btn.clicked.connect(self.canvas.reset_rotation)
        controls.addWidget(reset_btn)

        controls.addStretch()

        layout.addLayout(controls)

    def set_structure(self, structure: MoleculeStructure):
        """Display a molecule structure."""
        self.canvas.set_structure(structure)

    def clear(self):
        """Clear the viewer."""
        self.canvas._clear_plot()
        self.canvas.draw()
