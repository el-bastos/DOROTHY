"""
2D Molecule viewer widget using matplotlib.

Supports interactive atom picking for plane definition.
"""

import numpy as np
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
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
    """Canvas for displaying 2D molecule structure with atom picking."""

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

        # Connect pick event
        self.mpl_connect('button_press_event', self._on_click)

        self._clear_plot()

    def _clear_plot(self):
        """Clear the plot."""
        self.ax.clear()
        self.ax.set_aspect('equal')
        self.ax.axis('off')
        self.fig.tight_layout()

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

        # Project to 2D (use x, y coordinates - top-down view)
        # Could add rotation controls later
        x = coords[:, 0]
        y = coords[:, 1]

        # Center the molecule
        x = x - np.mean(x)
        y = y - np.mean(y)

        # Find bonds based on distance
        bonds = self._find_bonds(coords, symbols)

        # Draw bonds first (so they're behind atoms)
        for i, j in bonds:
            self.ax.plot([x[i], x[j]], [y[i], y[j]],
                        color='#404040', linewidth=2, zorder=1)

        # Store atom positions for hit testing
        self._atom_positions = list(zip(x, y))

        # Draw atoms
        for idx, (xi, yi, sym) in enumerate(zip(x, y, symbols)):
            color = ELEMENT_COLORS.get(sym, '#808080')
            # Size based on element (H smaller)
            size = 150 if sym == 'H' else 300

            # Highlight selected atoms
            if idx in self._selected_indices:
                edgecolor = '#FF0000'
                linewidth = 3
            else:
                edgecolor = '#404040'
                linewidth = 1

            # Draw atom circle
            self.ax.scatter(xi, yi, s=size, c=color, edgecolors=edgecolor,
                           linewidths=linewidth, zorder=2)

            # Label (skip H for cleaner view)
            if sym != 'H':
                self.ax.annotate(sym, (xi, yi), ha='center', va='center',
                               fontsize=8, fontweight='bold', zorder=3)

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
        """Handle mouse click for atom picking."""
        if not self._pick_enabled or event.inaxes != self.ax:
            return

        if not self._atom_positions:
            return

        # Only process left button clicks
        if event.button != 1:
            return

        # Find closest atom to click position using display coordinates
        click_x, click_y = event.x, event.y
        min_dist = float('inf')
        closest_idx = -1

        for idx, (ax, ay) in enumerate(self._atom_positions):
            # Transform data coords to display coords
            x_disp, y_disp = self.ax.transData.transform((ax, ay))
            dist = np.sqrt((x_disp - click_x) ** 2 + (y_disp - click_y) ** 2)
            if dist < min_dist:
                min_dist = dist
                closest_idx = idx

        # Check if click is close enough to an atom (threshold in pixels)
        if closest_idx >= 0 and min_dist < 25:  # 25 pixels threshold
            self.atom_picked.emit(closest_idx)

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


class MoleculeViewer(QWidget):
    """Widget for viewing molecule structures."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.canvas = MoleculeCanvas(self)
        layout.addWidget(self.canvas)

    def set_structure(self, structure: MoleculeStructure):
        """Display a molecule structure."""
        self.canvas.set_structure(structure)

    def clear(self):
        """Clear the viewer."""
        self.canvas._clear_plot()
        self.canvas.draw()
