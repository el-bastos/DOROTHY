"""
3D Interactive Slice Explorer widget using matplotlib.

Displays density slices stacked in 3D with interactive navigation.
Supports atom picking for plane definition, zoom controls, and bond visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel,
    QPushButton, QComboBox, QCheckBox, QSizePolicy, QGroupBox
)
from PyQt6.QtCore import Qt, pyqtSignal
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as mcolors

from dorothy.core.density import DensityCube

# Covalent radii (Å) for bond detection
COVALENT_RADII = {
    1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57,
    15: 1.07, 16: 1.05, 17: 1.02, 35: 1.20, 53: 1.39,
}


class SliceExplorerCanvas(FigureCanvas):
    """3D Canvas for displaying stacked density slices with atom picking."""

    # Signal emitted when an atom is clicked (passes atom index)
    atom_picked = pyqtSignal(int)

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

        # Atom selection state
        self._selected_indices: list[int] = []
        self._pick_enabled = False
        self._atom_positions_3d: list[tuple[float, float, float]] = []

        # Custom plane for slicing
        self._custom_plane_normal: np.ndarray | None = None

        # Zoom state
        self._zoom_level = 1.0
        self._base_limits = None  # Stores original axis limits

        # Show bonds when deformation density is available
        self._show_bonds = True

        # Connect mouse events
        self.mpl_connect('button_press_event', self._on_click)
        self.mpl_connect('scroll_event', self._on_scroll)

        self._setup_3d_axes()

    def _setup_3d_axes(self):
        """Set up the 3D axes."""
        self.fig.clear()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_facecolor('#f8f8f8')

        # Enable mouse interaction for rotation
        # This ensures drag-to-rotate works
        self.ax.mouse_init()

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

        # Get slices - use custom plane if defined
        basis_info = None
        if self._custom_plane_normal is not None:
            slices, basis_info = self._density_cube.get_oriented_slices(
                self._n_slices,
                self._custom_plane_normal
            )
            use_custom_plane = True
        else:
            slices = self._density_cube.get_z_slices(self._n_slices)
            use_custom_plane = False

        # Get grid extent in Angstrom
        origin = self._density_cube.origin * 0.529177
        nx, ny, nz = self._density_cube.shape
        axes = self._density_cube.axes * 0.529177

        x_extent = nx * axes[0, 0]
        y_extent = ny * axes[1, 1]

        # Determine contour levels based on density type
        if self._density_type == "deformation":
            all_data = np.concatenate([s.flatten() for _, s in slices])
            max_abs = max(abs(all_data.min()), abs(all_data.max()))
            if max_abs > 0:
                levels_pos = [0.005, 0.01, 0.02, 0.04, 0.08]
                levels_neg = [-0.08, -0.04, -0.02, -0.01, -0.005]
            else:
                levels_pos = []
                levels_neg = []
        else:
            all_data = np.concatenate([s.flatten() for _, s in slices])
            max_val = all_data.max()
            if max_val > 0:
                levels_pos = np.linspace(max_val * 0.05, max_val * 0.8, 5)
            else:
                levels_pos = []
            levels_neg = []

        # Draw each slice
        if use_custom_plane and basis_info:
            self._draw_oriented_slices(slices, basis_info, levels_pos, levels_neg)
            # Calculate limits for oriented slices
            center = basis_info['center']
            half = basis_info['half_extent'] + 1.0
            xlim = (center[0] - half, center[0] + half)
            ylim = (center[1] - half, center[1] + half)
            zlim = (center[2] - half, center[2] + half)
        else:
            self._draw_z_slices(slices, origin, x_extent, y_extent, levels_pos, levels_neg)
            # Calculate limits for Z-axis slices
            x = np.linspace(origin[0], origin[0] + x_extent, slices[0][1].shape[0])
            y = np.linspace(origin[1], origin[1] + y_extent, slices[0][1].shape[1])
            z_coords = [z for z, _ in slices]
            padding = 0.5
            xlim = (x.min() - padding, x.max() + padding)
            ylim = (y.min() - padding, y.max() + padding)
            zlim = (min(z_coords) - padding, max(z_coords) + padding)

        # Draw bonds and atoms
        self._draw_bonds()
        self._draw_atoms()

        # Store base limits for zoom
        self._base_limits = (xlim, ylim, zlim)

        # Apply current zoom level
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.ax.set_zlim(zlim)
        self._apply_zoom()

        self.fig.tight_layout()
        self.draw()

    def _draw_z_slices(self, slices, origin, x_extent, y_extent, levels_pos, levels_neg):
        """Draw slices along Z-axis (original behavior)."""
        x = np.linspace(origin[0], origin[0] + x_extent, slices[0][1].shape[0])
        y = np.linspace(origin[1], origin[1] + y_extent, slices[0][1].shape[1])
        X, Y = np.meshgrid(x, y, indexing='ij')

        for i, (z_coord, slice_data) in enumerate(slices):
            is_highlighted = (i == self._highlighted_slice)
            alpha = 0.9 if is_highlighted else 0.15
            Z = np.full_like(X, z_coord)

            if self._show_contours and len(levels_pos) > 0:
                self._draw_slice_contours(
                    X, Y, Z, slice_data, z_coord,
                    levels_pos, levels_neg,
                    alpha, is_highlighted
                )
            else:
                self._draw_slice_plane(X, Y, z_coord, slice_data, alpha, is_highlighted)

    def _draw_oriented_slices(self, slices, basis_info, levels_pos, levels_neg):
        """Draw slices along custom plane orientation."""
        normal = basis_info['normal']
        u = basis_info['u']
        v = basis_info['v']
        center = basis_info['center']
        half_extent = basis_info['half_extent']

        slice_shape = slices[0][1].shape
        u_range = np.linspace(-half_extent, half_extent, slice_shape[0])
        v_range = np.linspace(-half_extent, half_extent, slice_shape[1])
        U, V = np.meshgrid(u_range, v_range, indexing='ij')

        for i, (offset, slice_data) in enumerate(slices):
            is_highlighted = (i == self._highlighted_slice)
            alpha = 0.9 if is_highlighted else 0.15

            # Calculate 3D position of each point in the slice
            slice_center = center + normal * offset
            X = slice_center[0] + U * u[0] + V * v[0]
            Y = slice_center[1] + U * u[1] + V * v[1]
            Z = slice_center[2] + U * u[2] + V * v[2]

            if self._show_contours and len(levels_pos) > 0:
                self._draw_oriented_contours(
                    X, Y, Z, slice_data,
                    levels_pos, levels_neg,
                    alpha, is_highlighted
                )
            else:
                self._draw_oriented_plane(X, Y, Z, alpha, is_highlighted)

    def _draw_oriented_contours(self, X, Y, Z, slice_data, levels_pos, levels_neg, alpha, is_highlighted):
        """Draw contour lines on an oriented slice."""
        # Draw the plane first
        self._draw_oriented_plane(X, Y, Z, alpha * 0.3, is_highlighted)

        # Draw positive contours using surface plot approach
        if len(levels_pos) > 0 and slice_data.max() > min(levels_pos):
            try:
                # For oriented slices, we can't use offset parameter
                # Instead, we plot the contours in 3D directly
                import matplotlib.pyplot as plt
                from matplotlib import cm

                # Get contour lines in 2D, then map to 3D
                fig_temp = plt.figure()
                ax_temp = fig_temp.add_subplot(111)
                cs = ax_temp.contour(slice_data, levels=levels_pos)
                plt.close(fig_temp)

                for collection in cs.collections:
                    for path in collection.get_paths():
                        vertices = path.vertices
                        if len(vertices) > 1:
                            # Map 2D contour indices to 3D coordinates
                            ni, nj = slice_data.shape
                            xs = []
                            ys = []
                            zs = []
                            for vi, vj in vertices:
                                # vi, vj are in array index space
                                ii = int(np.clip(vi, 0, ni - 1))
                                jj = int(np.clip(vj, 0, nj - 1))
                                xs.append(X[ii, jj])
                                ys.append(Y[ii, jj])
                                zs.append(Z[ii, jj])
                            self.ax.plot(xs, ys, zs,
                                       color='#000000' if self._color_mode == 'bw' else '#2563eb',
                                       linewidth=2.0 if is_highlighted else 0.6,
                                       alpha=1.0 if is_highlighted else 0.4)
            except Exception:
                pass

        # Draw negative contours
        if len(levels_neg) > 0 and slice_data.min() < max(levels_neg):
            try:
                fig_temp = plt.figure()
                ax_temp = fig_temp.add_subplot(111)
                cs = ax_temp.contour(slice_data, levels=levels_neg)
                plt.close(fig_temp)

                for collection in cs.collections:
                    for path in collection.get_paths():
                        vertices = path.vertices
                        if len(vertices) > 1:
                            ni, nj = slice_data.shape
                            xs = []
                            ys = []
                            zs = []
                            for vi, vj in vertices:
                                ii = int(np.clip(vi, 0, ni - 1))
                                jj = int(np.clip(vj, 0, nj - 1))
                                xs.append(X[ii, jj])
                                ys.append(Y[ii, jj])
                                zs.append(Z[ii, jj])
                            self.ax.plot(xs, ys, zs,
                                       color='#666666' if self._color_mode == 'bw' else '#dc2626',
                                       linestyle='--',
                                       linewidth=2.0 if is_highlighted else 0.6,
                                       alpha=1.0 if is_highlighted else 0.4)
            except Exception:
                pass

    def _draw_oriented_plane(self, X, Y, Z, alpha, is_highlighted):
        """Draw a semi-transparent oriented slice plane."""
        # Get corners of the slice
        corners = [
            [X[0, 0], Y[0, 0], Z[0, 0]],
            [X[-1, 0], Y[-1, 0], Z[-1, 0]],
            [X[-1, -1], Y[-1, -1], Z[-1, -1]],
            [X[0, -1], Y[0, -1], Z[0, -1]],
        ]

        if is_highlighted:
            color = '#2563eb' if self._color_mode == 'color' else '#666666'
            face_alpha = 0.3
        else:
            color = '#cccccc'
            face_alpha = 0.1

        poly = Poly3DCollection([corners], alpha=face_alpha * alpha,
                                facecolor=color, edgecolor='#999999', linewidth=0.5)
        self.ax.add_collection3d(poly)

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

    def _draw_bonds(self):
        """Draw bonds between atoms in 3D."""
        if not self._density_cube or not self._density_cube.atoms:
            return

        if not self._show_bonds:
            return

        atoms = self._density_cube.atoms
        n_atoms = len(atoms)

        for i in range(n_atoms):
            z1, x1, y1, z1_coord = atoms[i]
            x1 *= 0.529177
            y1 *= 0.529177
            z1_coord *= 0.529177
            r1 = COVALENT_RADII.get(z1, 1.5)

            for j in range(i + 1, n_atoms):
                z2, x2, y2, z2_coord = atoms[j]
                x2 *= 0.529177
                y2 *= 0.529177
                z2_coord *= 0.529177
                r2 = COVALENT_RADII.get(z2, 1.5)

                dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2_coord - z1_coord)**2)
                max_bond_dist = (r1 + r2) * 1.3

                if dist < max_bond_dist:
                    self.ax.plot(
                        [x1, x2], [y1, y2], [z1_coord, z2_coord],
                        color='#404040',
                        linewidth=2,
                        alpha=0.7,
                        zorder=1
                    )

    def _draw_atoms(self):
        """Draw atom positions in 3D with selection highlighting."""
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

        # Store positions for hit testing
        self._atom_positions_3d = []

        for idx, atom in enumerate(self._density_cube.atoms):
            z_num, ax, ay, az = atom
            # Convert from Bohr to Angstrom
            x = ax * 0.529177
            y = ay * 0.529177
            z = az * 0.529177

            self._atom_positions_3d.append((x, y, z))

            color = element_colors.get(z_num, '#808080')
            size = 20 if z_num == 1 else 50

            # Highlight selected atoms
            if idx in self._selected_indices:
                edgecolor = '#FF0000'
                linewidth = 2
                size *= 1.5
            else:
                edgecolor = '#404040'
                linewidth = 0.5

            self.ax.scatter([x], [y], [z],
                          c=color, s=size,
                          edgecolors=edgecolor,
                          linewidths=linewidth,
                          alpha=0.9,
                          depthshade=True,
                          zorder=10)

    def _on_click(self, event):
        """Handle mouse click for atom picking in 3D.

        Note: We use 'button_press_event' but only handle atom picking.
        Matplotlib's default 3D rotation (drag) still works because we
        don't block the event - we just emit a signal if an atom is clicked.
        """
        # Only process left button clicks for picking
        if event.button != 1:  # 1 = left button
            return

        if not self._pick_enabled or event.inaxes != self.ax:
            return

        if not self._atom_positions_3d:
            return

        # Get click position in display coordinates
        click_x, click_y = event.x, event.y
        if click_x is None or click_y is None:
            return

        # Project all atoms to display coordinates and find closest
        from mpl_toolkits.mplot3d import proj3d

        min_dist = float('inf')
        closest_idx = -1

        for idx, (ax, ay, az) in enumerate(self._atom_positions_3d):
            # Project 3D point to display coordinates
            x2d, y2d, _ = proj3d.proj_transform(ax, ay, az, self.ax.get_proj())
            # Convert to display coordinates
            x_disp, y_disp = self.ax.transData.transform((x2d, y2d))
            dist = np.sqrt((x_disp - click_x) ** 2 + (y_disp - click_y) ** 2)
            if dist < min_dist:
                min_dist = dist
                closest_idx = idx

        # Threshold for click detection (in display/pixel units)
        if closest_idx >= 0 and min_dist < 30:  # 30 pixels threshold
            self.atom_picked.emit(closest_idx)

    def _on_scroll(self, event):
        """Handle scroll wheel for zoom."""
        if event.inaxes != self.ax:
            return

        if event.button == 'up':
            self._zoom_level *= 0.9  # Zoom in
        else:
            self._zoom_level *= 1.1  # Zoom out

        self._zoom_level = max(0.2, min(5.0, self._zoom_level))
        self._apply_zoom()

    def _apply_zoom(self):
        """Apply current zoom level to axis limits."""
        if self._base_limits is None:
            return

        xlim, ylim, zlim = self._base_limits

        # Calculate center
        cx = (xlim[0] + xlim[1]) / 2
        cy = (ylim[0] + ylim[1]) / 2
        cz = (zlim[0] + zlim[1]) / 2

        # Apply zoom
        half_x = (xlim[1] - xlim[0]) / 2 * self._zoom_level
        half_y = (ylim[1] - ylim[0]) / 2 * self._zoom_level
        half_z = (zlim[1] - zlim[0]) / 2 * self._zoom_level

        self.ax.set_xlim(cx - half_x, cx + half_x)
        self.ax.set_ylim(cy - half_y, cy + half_y)
        self.ax.set_zlim(cz - half_z, cz + half_z)

        self.draw()

    def zoom_in(self):
        """Zoom in by 10%."""
        self._zoom_level *= 0.9
        self._zoom_level = max(0.2, self._zoom_level)
        self._apply_zoom()

    def zoom_out(self):
        """Zoom out by 10%."""
        self._zoom_level *= 1.1
        self._zoom_level = min(5.0, self._zoom_level)
        self._apply_zoom()

    def reset_zoom(self):
        """Reset zoom to default."""
        self._zoom_level = 1.0
        self._apply_zoom()

    def rotate_left(self):
        """Rotate view left by 15 degrees."""
        elev, azim = self.ax.elev, self.ax.azim
        self.ax.view_init(elev=elev, azim=azim - 15)
        self.draw()

    def rotate_right(self):
        """Rotate view right by 15 degrees."""
        elev, azim = self.ax.elev, self.ax.azim
        self.ax.view_init(elev=elev, azim=azim + 15)
        self.draw()

    def rotate_up(self):
        """Rotate view up by 15 degrees."""
        elev, azim = self.ax.elev, self.ax.azim
        new_elev = min(90, elev + 15)
        self.ax.view_init(elev=new_elev, azim=azim)
        self.draw()

    def rotate_down(self):
        """Rotate view down by 15 degrees."""
        elev, azim = self.ax.elev, self.ax.azim
        new_elev = max(-90, elev - 15)
        self.ax.view_init(elev=new_elev, azim=azim)
        self.draw()

    def set_selection(self, indices: list[int]):
        """Update visual selection."""
        self._selected_indices = indices.copy()
        if self._density_cube:
            self._draw_slices()

    def set_pick_enabled(self, enabled: bool):
        """Enable or disable atom picking."""
        self._pick_enabled = enabled

    def set_custom_plane(self, normal: np.ndarray | None):
        """Set custom slicing plane and redraw."""
        self._custom_plane_normal = normal
        if self._density_cube:
            self._draw_slices()

    def set_show_bonds(self, show: bool):
        """Enable or disable bond visualization."""
        self._show_bonds = show
        if self._density_cube:
            self._draw_slices()


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

        self.bonds_check = QCheckBox("Bonds")
        self.bonds_check.setChecked(True)
        self.bonds_check.stateChanged.connect(self._on_bonds_toggled)
        contour_group.addWidget(self.bonds_check)
        controls.addLayout(contour_group)

        # Zoom controls
        zoom_group = QVBoxLayout()
        zoom_label = QLabel("Zoom:")
        zoom_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        zoom_group.addWidget(zoom_label)

        zoom_row = QHBoxLayout()
        zoom_in_btn = QPushButton("+")
        zoom_in_btn.setMaximumWidth(30)
        zoom_in_btn.clicked.connect(self.canvas.zoom_in)
        zoom_row.addWidget(zoom_in_btn)

        zoom_out_btn = QPushButton("-")
        zoom_out_btn.setMaximumWidth(30)
        zoom_out_btn.clicked.connect(self.canvas.zoom_out)
        zoom_row.addWidget(zoom_out_btn)

        zoom_reset_btn = QPushButton("Reset")
        zoom_reset_btn.clicked.connect(self.canvas.reset_zoom)
        zoom_row.addWidget(zoom_reset_btn)

        zoom_group.addLayout(zoom_row)
        controls.addLayout(zoom_group)

        # Rotation controls
        rotate_group = QVBoxLayout()
        rotate_label = QLabel("Rotate:")
        rotate_label.setStyleSheet("font-weight: bold; font-size: 10pt;")
        rotate_group.addWidget(rotate_label)

        rotate_row = QHBoxLayout()
        rotate_left_btn = QPushButton("\u2190")  # Left arrow
        rotate_left_btn.setMaximumWidth(30)
        rotate_left_btn.setToolTip("Rotate left (or drag mouse)")
        rotate_left_btn.clicked.connect(self.canvas.rotate_left)
        rotate_row.addWidget(rotate_left_btn)

        rotate_right_btn = QPushButton("\u2192")  # Right arrow
        rotate_right_btn.setMaximumWidth(30)
        rotate_right_btn.setToolTip("Rotate right (or drag mouse)")
        rotate_right_btn.clicked.connect(self.canvas.rotate_right)
        rotate_row.addWidget(rotate_right_btn)

        rotate_up_btn = QPushButton("\u2191")  # Up arrow
        rotate_up_btn.setMaximumWidth(30)
        rotate_up_btn.setToolTip("Rotate up")
        rotate_up_btn.clicked.connect(self.canvas.rotate_up)
        rotate_row.addWidget(rotate_up_btn)

        rotate_down_btn = QPushButton("\u2193")  # Down arrow
        rotate_down_btn.setMaximumWidth(30)
        rotate_down_btn.setToolTip("Rotate down")
        rotate_down_btn.clicked.connect(self.canvas.rotate_down)
        rotate_row.addWidget(rotate_down_btn)

        rotate_group.addLayout(rotate_row)
        controls.addLayout(rotate_group)

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

    def _on_bonds_toggled(self, state: int):
        """Handle bonds checkbox toggle."""
        self.canvas.set_show_bonds(state == Qt.CheckState.Checked.value)

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
