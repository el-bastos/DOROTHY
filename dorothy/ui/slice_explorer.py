"""
3D Interactive Slice Explorer widget using matplotlib.

Displays density slices stacked in 3D with interactive navigation.
"""

import numpy as np
import matplotlib.pyplot as plt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel,
    QPushButton, QComboBox, QCheckBox, QSizePolicy, QGroupBox
)
from PyQt6.QtCore import Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as mcolors

from dorothy.core.density import DensityCube


class SliceExplorerCanvas(FigureCanvas):
    """3D Canvas for displaying stacked density slices."""

    def __init__(self, parent=None):
        self.fig = Figure(figsize=(6, 6), dpi=100)
        self.fig.patch.set_facecolor('#f8f8f8')
        super().__init__(self.fig)
        self.setParent(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._density_cube: DensityCube | None = None
        self._n_slices = 15
        self._highlighted_slice = 7  # Middle slice
        self._show_contours = True
        self._density_type = "promolecule"  # or "deformation"
        self._color_mode = "bw"  # or "color"

        self._setup_3d_axes()

    def _setup_3d_axes(self):
        """Set up the 3D axes."""
        self.fig.clear()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_facecolor('#f8f8f8')

        # Set initial view angle
        self.ax.view_init(elev=25, azim=-60)

        # Clean axis appearance
        self.ax.set_xlabel('X (Å)', fontsize=9)
        self.ax.set_ylabel('Y (Å)', fontsize=9)
        self.ax.set_zlabel('Z (Å)', fontsize=9)

        self.fig.tight_layout()

    def set_density_cube(self, cube: DensityCube, density_type: str = "promolecule"):
        """Set the density cube to visualize."""
        self._density_cube = cube
        self._density_type = density_type
        self._highlighted_slice = self._n_slices // 2
        self._draw_slices()

    def set_n_slices(self, n: int):
        """Set number of slices to display."""
        self._n_slices = n
        self._highlighted_slice = min(self._highlighted_slice, n - 1)
        if self._density_cube:
            self._draw_slices()

    def set_highlighted_slice(self, index: int):
        """Set which slice is highlighted."""
        self._highlighted_slice = max(0, min(index, self._n_slices - 1))
        if self._density_cube:
            self._draw_slices()

    def set_color_mode(self, mode: str):
        """Set color mode ('bw' or 'color')."""
        self._color_mode = mode
        if self._density_cube:
            self._draw_slices()

    def set_show_contours(self, show: bool):
        """Toggle contour display."""
        self._show_contours = show
        if self._density_cube:
            self._draw_slices()

    def _draw_slices(self):
        """Draw all slices in 3D."""
        self._setup_3d_axes()

        if not self._density_cube:
            self.ax.text2D(0.5, 0.5, 'No density data',
                          ha='center', va='center', transform=self.ax.transAxes,
                          fontsize=12, color='#888')
            self.draw()
            return

        # Get slices
        slices = self._density_cube.get_z_slices(self._n_slices)

        # Get grid extent in Angstrom
        origin = self._density_cube.origin * 0.529177
        nx, ny, nz = self._density_cube.shape
        axes = self._density_cube.axes * 0.529177

        x_extent = nx * axes[0, 0]
        y_extent = ny * axes[1, 1]

        # Create x, y coordinates for slices
        x = np.linspace(origin[0], origin[0] + x_extent, slices[0][1].shape[0])
        y = np.linspace(origin[1], origin[1] + y_extent, slices[0][1].shape[1])
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Determine contour levels based on density type
        if self._density_type == "deformation":
            # Deformation density has positive and negative values
            all_data = np.concatenate([s.flatten() for _, s in slices])
            max_abs = max(abs(all_data.min()), abs(all_data.max()))
            if max_abs > 0:
                levels_pos = [0.005, 0.01, 0.02, 0.04, 0.08]
                levels_neg = [-0.08, -0.04, -0.02, -0.01, -0.005]
            else:
                levels_pos = []
                levels_neg = []
        else:
            # Promolecule density is always positive
            all_data = np.concatenate([s.flatten() for _, s in slices])
            max_val = all_data.max()
            if max_val > 0:
                # Use fixed levels that make sense for electron density
                levels_pos = np.linspace(max_val * 0.05, max_val * 0.8, 5)
            else:
                levels_pos = []
            levels_neg = []

        # Draw each slice
        for i, (z_coord, slice_data) in enumerate(slices):
            is_highlighted = (i == self._highlighted_slice)
            alpha = 0.9 if is_highlighted else 0.15

            # Create the slice plane
            Z = np.full_like(X, z_coord)

            if self._show_contours and len(levels_pos) > 0:
                # Draw contours on the slice
                self._draw_slice_contours(
                    X, Y, Z, slice_data, z_coord,
                    levels_pos, levels_neg,
                    alpha, is_highlighted
                )
            else:
                # Just draw a semi-transparent plane
                self._draw_slice_plane(X, Y, z_coord, slice_data, alpha, is_highlighted)

        # Draw atoms if available
        self._draw_atoms()

        # Set axis limits
        z_coords = [z for z, _ in slices]
        padding = 0.5
        self.ax.set_xlim(x.min() - padding, x.max() + padding)
        self.ax.set_ylim(y.min() - padding, y.max() + padding)
        self.ax.set_zlim(min(z_coords) - padding, max(z_coords) + padding)

        self.fig.tight_layout()
        self.draw()

    def _draw_slice_contours(self, X, Y, Z, slice_data, z_coord,
                             levels_pos, levels_neg, alpha, is_highlighted):
        """Draw contour lines on a slice at given z-coordinate."""
        # First draw the plane behind the contours
        self._draw_slice_plane(X, Y, z_coord, slice_data, alpha * 0.3, is_highlighted)

        # For 3D contour with zdir='z', we pass X, Y grids and data directly
        # The data shape should match (nx, ny) where X is (nx, ny) and Y is (nx, ny)

        # Draw positive contours
        if len(levels_pos) > 0 and slice_data.max() > min(levels_pos):
            try:
                self.ax.contour(
                    X, Y, slice_data,
                    levels=levels_pos,
                    zdir='z',
                    offset=z_coord,
                    colors='#000000' if self._color_mode == 'bw' else '#2563eb',
                    linewidths=2.0 if is_highlighted else 0.6,
                    alpha=1.0 if is_highlighted else 0.4
                )
            except Exception:
                pass

        # Draw negative contours (for deformation density)
        if len(levels_neg) > 0 and slice_data.min() < max(levels_neg):
            try:
                self.ax.contour(
                    X, Y, slice_data,
                    levels=levels_neg,
                    zdir='z',
                    offset=z_coord,
                    colors='#666666' if self._color_mode == 'bw' else '#dc2626',
                    linestyles='dashed',
                    linewidths=2.0 if is_highlighted else 0.6,
                    alpha=1.0 if is_highlighted else 0.4
                )
            except Exception:
                pass

    def _draw_slice_plane(self, X, Y, z_coord, slice_data, alpha, is_highlighted):
        """Draw a semi-transparent slice plane as a simple rectangle."""
        # Get bounds
        x_min, x_max = X.min(), X.max()
        y_min, y_max = Y.min(), Y.max()

        # Draw a simple colored rectangle for the slice plane
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        # Define the corners of the rectangle
        verts = [
            [(x_min, y_min, z_coord),
             (x_max, y_min, z_coord),
             (x_max, y_max, z_coord),
             (x_min, y_max, z_coord)]
        ]

        # Choose color based on density
        if is_highlighted:
            color = '#2563eb' if self._color_mode == 'color' else '#666666'
            face_alpha = 0.3
        else:
            color = '#cccccc'
            face_alpha = 0.1

        poly = Poly3DCollection(verts, alpha=face_alpha, facecolor=color,
                                edgecolor='#999999', linewidth=0.5)
        self.ax.add_collection3d(poly)

    def _draw_atoms(self):
        """Draw atom positions in 3D."""
        if not self._density_cube or not self._density_cube.atoms:
            return

        # Element colors
        element_colors = {
            1: '#FFFFFF',   # H
            6: '#909090',   # C
            7: '#3050F8',   # N
            8: '#FF0D0D',   # O
            9: '#90E050',   # F
            15: '#FF8000',  # P
            16: '#FFFF30',  # S
            17: '#1FF01F',  # Cl
            35: '#A62929',  # Br
            53: '#940094',  # I
        }

        for atom in self._density_cube.atoms:
            z_num, ax, ay, az = atom
            # Convert from Bohr to Angstrom
            x = ax * 0.529177
            y = ay * 0.529177
            z = az * 0.529177

            color = element_colors.get(z_num, '#808080')
            size = 20 if z_num == 1 else 50

            self.ax.scatter([x], [y], [z],
                          c=color, s=size,
                          edgecolors='#404040',
                          linewidths=0.5,
                          alpha=0.8,
                          depthshade=True)


class SliceExplorer(QWidget):
    """Interactive 3D Slice Explorer widget."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._promolecule_cube: DensityCube | None = None
        self._deformation_cube: DensityCube | None = None

        self._setup_ui()

    def _setup_ui(self):
        """Set up the widget UI."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)

        # 3D Canvas
        self.canvas = SliceExplorerCanvas(self)
        self.canvas.setMinimumSize(400, 400)
        layout.addWidget(self.canvas, stretch=1)

        # Controls
        controls = QHBoxLayout()
        controls.setSpacing(15)

        # Density type selector
        type_group = QVBoxLayout()
        type_label = QLabel("Density:")
        type_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        type_group.addWidget(type_label)

        self.density_combo = QComboBox()
        self.density_combo.addItems(["Promolecule", "Deformation"])
        self.density_combo.currentIndexChanged.connect(self._on_density_type_changed)
        type_group.addWidget(self.density_combo)
        controls.addLayout(type_group)

        # Slice slider
        slider_group = QVBoxLayout()
        slider_label = QLabel("Slice:")
        slider_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        slider_group.addWidget(slider_label)

        slider_row = QHBoxLayout()
        self.slice_slider = QSlider(Qt.Orientation.Horizontal)
        self.slice_slider.setMinimum(0)
        self.slice_slider.setMaximum(14)
        self.slice_slider.setValue(7)
        self.slice_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.slice_slider.setTickInterval(1)
        self.slice_slider.valueChanged.connect(self._on_slice_changed)
        slider_row.addWidget(self.slice_slider, stretch=1)

        self.slice_label = QLabel("8/15")
        self.slice_label.setMinimumWidth(45)
        slider_row.addWidget(self.slice_label)
        slider_group.addLayout(slider_row)
        controls.addLayout(slider_group, stretch=1)

        # Color mode
        color_group = QVBoxLayout()
        color_label = QLabel("Color:")
        color_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        color_group.addWidget(color_label)

        self.color_combo = QComboBox()
        self.color_combo.addItems(["B&W", "Color"])
        self.color_combo.currentIndexChanged.connect(self._on_color_mode_changed)
        color_group.addWidget(self.color_combo)
        controls.addLayout(color_group)

        # Contours toggle
        contour_group = QVBoxLayout()
        contour_label = QLabel("Display:")
        contour_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        contour_group.addWidget(contour_label)

        self.contours_check = QCheckBox("Contours")
        self.contours_check.setChecked(True)
        self.contours_check.stateChanged.connect(self._on_contours_toggled)
        contour_group.addWidget(self.contours_check)
        controls.addLayout(contour_group)

        layout.addLayout(controls)

    def set_density_cubes(self, promolecule: DensityCube | None = None,
                          deformation: DensityCube | None = None,
                          n_slices: int = 15):
        """Set the density cubes to visualize."""
        self._promolecule_cube = promolecule
        self._deformation_cube = deformation

        # Update slice slider
        self.slice_slider.setMaximum(n_slices - 1)
        self.slice_slider.setValue(n_slices // 2)
        self.canvas.set_n_slices(n_slices)

        # Update density type combo
        self.density_combo.clear()
        if promolecule is not None:
            self.density_combo.addItem("Promolecule")
        if deformation is not None:
            self.density_combo.addItem("Deformation")

        # Show first available density
        if promolecule is not None:
            self.canvas.set_density_cube(promolecule, "promolecule")
        elif deformation is not None:
            self.canvas.set_density_cube(deformation, "deformation")

        self._update_slice_label()

    def _on_density_type_changed(self, index: int):
        """Handle density type selection change."""
        text = self.density_combo.currentText().lower()
        if text == "promolecule" and self._promolecule_cube:
            self.canvas.set_density_cube(self._promolecule_cube, "promolecule")
        elif text == "deformation" and self._deformation_cube:
            self.canvas.set_density_cube(self._deformation_cube, "deformation")

    def _on_slice_changed(self, value: int):
        """Handle slice slider change."""
        self.canvas.set_highlighted_slice(value)
        self._update_slice_label()

    def _on_color_mode_changed(self, index: int):
        """Handle color mode change."""
        mode = "bw" if index == 0 else "color"
        self.canvas.set_color_mode(mode)

    def _on_contours_toggled(self, state: int):
        """Handle contours checkbox toggle."""
        self.canvas.set_show_contours(state == Qt.CheckState.Checked.value)

    def _update_slice_label(self):
        """Update the slice number label."""
        current = self.slice_slider.value() + 1
        total = self.slice_slider.maximum() + 1
        self.slice_label.setText(f"{current}/{total}")

    def clear(self):
        """Clear the explorer."""
        self._promolecule_cube = None
        self._deformation_cube = None
        self.canvas._density_cube = None
        self.canvas._setup_3d_axes()
        self.canvas.draw()
