"""
3D Interactive Slice Explorer widget using matplotlib.

Displays density slices stacked in 3D with interactive navigation.
Supports atom picking for plane definition, zoom controls, and bond visualization.
"""

import numpy as np
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel,
    QPushButton, QComboBox, QCheckBox, QSizePolicy, QButtonGroup,
    QFrame, QDialog, QDialogButtonBox
)
from PyQt6.QtCore import Qt, pyqtSignal, QTimer
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.colors import Normalize
from matplotlib import cm

from dorothy.core.density import DensityCube
from dorothy.core.constants import (
    BOHR_TO_ANGSTROM,
    Z_TO_SYMBOL,
    COVALENT_RADII,
    ELEMENT_COLORS,
    CONTOUR_LEVELS_BASE,
    SLICE_ATOM_THICKNESS,
)


class Slice2DPreviewDialog(QDialog):
    """Dialog showing a 2D contour plot of the current slice."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("2D Slice Preview")
        self.setMinimumSize(500, 500)

        layout = QVBoxLayout(self)

        # Matplotlib figure
        self.fig = Figure(figsize=(6, 6), dpi=100)
        self.fig.patch.set_facecolor('white')
        self.canvas = FigureCanvas(self.fig)
        layout.addWidget(self.canvas)

        # Info label
        self.info_label = QLabel()
        self.info_label.setStyleSheet("font-size: 10pt; color: #666;")
        layout.addWidget(self.info_label)

        # Close button
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)

    def show_slice(self, slice_data: np.ndarray, z_coord: float, density_type: str,
                   color_mode: str, atoms: list = None, contour_scale: float = 1.0,
                   origin: np.ndarray = None, axes: np.ndarray = None,
                   render_mode: str = "contour", colormap: str = "RdYlBu_r"):
        """Display a 2D contour plot or heatmap of the slice.

        Args:
            slice_data: 2D array of density values
            z_coord: Z-coordinate of the slice in Angstrom
            density_type: "promolecule" or "deformation"
            color_mode: "bw" or "color"
            atoms: List of (Z, x, y, z) in Bohr - atoms near this slice will be shown
            contour_scale: Contour level scale factor
            origin: Grid origin in Bohr
            axes: Grid axes in Bohr
            render_mode: "contour" or "heatmap"
            colormap: Matplotlib colormap name for heatmap mode
        """
        self.fig.clear()
        ax = self.fig.add_subplot(111)

        # Create coordinate arrays if origin/axes provided
        if origin is not None and axes is not None:
            nx, ny = slice_data.shape
            x = np.linspace(origin[0], origin[0] + nx * axes[0, 0], nx) * BOHR_TO_ANGSTROM
            y = np.linspace(origin[1], origin[1] + ny * axes[1, 1], ny) * BOHR_TO_ANGSTROM
            X, Y = np.meshgrid(x, y, indexing='ij')
            extent = [x.min(), x.max(), y.min(), y.max()]
        else:
            X, Y = np.meshgrid(np.arange(slice_data.shape[0]),
                               np.arange(slice_data.shape[1]), indexing='ij')
            extent = None

        # Heatmap mode - use imshow
        if render_mode == "heatmap":
            if density_type == "deformation":
                max_abs = max(abs(slice_data.min()), abs(slice_data.max()))
                vmin, vmax = -max_abs, max_abs
                info_text = f"Deformation density at z = {z_coord:.2f} Å\nRed = accumulation | Blue = depletion"
            else:
                vmin, vmax = slice_data.min(), slice_data.max()
                info_text = f"Promolecule density at z = {z_coord:.2f} Å"

            # Note: imshow expects data in (rows, cols) = (y, x) format
            # Our data is (x, y), so we transpose
            im = ax.imshow(
                slice_data.T,
                cmap=colormap,
                vmin=vmin,
                vmax=vmax,
                extent=extent,
                origin='lower',
                aspect='equal',
                interpolation='bilinear'
            )
            self.fig.colorbar(im, ax=ax, shrink=0.8)
            self.info_label.setText(info_text)
        else:
            # Contour mode (original behavior)
            scale = contour_scale
            if density_type == "deformation":
                base_levels = np.array(CONTOUR_LEVELS_BASE)
                levels_pos = list(base_levels * scale)
                levels_neg = list(-base_levels[::-1] * scale)

                # Draw positive contours
                valid_pos = [l for l in levels_pos if l < slice_data.max()]
                if valid_pos:
                    color = 'black' if color_mode == 'bw' else '#2563eb'
                    ax.contour(X, Y, slice_data, levels=valid_pos, colors=color, linewidths=1.0)

                # Draw negative contours
                valid_neg = [l for l in levels_neg if l > slice_data.min()]
                if valid_neg:
                    color = '#666666' if color_mode == 'bw' else '#dc2626'
                    ax.contour(X, Y, slice_data, levels=valid_neg, colors=color,
                              linewidths=1.0, linestyles='dashed')

                self.info_label.setText(
                    f"Deformation density at z = {z_coord:.2f} Å\n"
                    "Solid = electron accumulation (bonds) | Dashed = depletion"
                )
            else:
                # Promolecule
                max_val = slice_data.max()
                if max_val > 0.001:
                    n_contours = max(3, int(7 / scale))
                    start = max_val * 0.02 * scale
                    end = max_val * 0.8
                    if start < end:
                        levels = list(np.linspace(start, end, n_contours))
                        color = 'black' if color_mode == 'bw' else '#2563eb'
                        ax.contour(X, Y, slice_data, levels=levels, colors=color, linewidths=1.0)

                self.info_label.setText(f"Promolecule density at z = {z_coord:.2f} Å")

        # Draw atoms near this slice
        if atoms and origin is not None:
            for atom in atoms:
                z_num, ax_b, ay_b, az_b = atom
                az = az_b * BOHR_TO_ANGSTROM
                if abs(az - z_coord) < SLICE_ATOM_THICKNESS:
                    ax_a = ax_b * BOHR_TO_ANGSTROM
                    ay_a = ay_b * BOHR_TO_ANGSTROM
                    color = ELEMENT_COLORS.get(z_num, '#808080')
                    size = 100 if z_num == 1 else 200
                    ax.scatter([ax_a], [ay_a], c=color, s=size, edgecolors='black',
                              linewidths=1, alpha=0.8, zorder=10)
                    # Label
                    symbol = Z_TO_SYMBOL.get(z_num, '?')
                    ax.annotate(symbol, (ax_a, ay_a), fontsize=8, ha='center', va='center')

        ax.set_aspect('equal')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_title(f'Slice at z = {z_coord:.2f} Å')

        self.fig.tight_layout()
        self.canvas.draw()


class SliceExplorerCanvas(FigureCanvas):
    """3D Canvas for displaying stacked density slices with atom picking."""

    # Signal emitted when an atom is clicked (passes atom index)
    atom_picked = pyqtSignal(int)
    # Signal emitted when a bond is clicked (passes tuple of atom indices and bond info)
    bond_clicked = pyqtSignal(int, int, float)  # atom1_idx, atom2_idx, bond_length

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

        # Zoom state
        self._zoom_level = 1.0
        self._base_limits = None  # Stores original axis limits

        # Show bonds when deformation density is available
        self._show_bonds = True

        # Show only highlighted slice (hide others)
        self._show_only_highlighted = False

        # Density hover display
        self._show_density_on_hover = False
        self._density_text = None

        # Contour level scale (1.0 = default, <1 = tighter/more, >1 = looser/fewer)
        self._contour_scale = 1.0

        # Render mode: "contour" or "heatmap"
        self._render_mode = "contour"

        # Colormap for heatmap mode
        self._heatmap_colormap = "RdYlBu_r"  # Red (high) -> Blue (low)

        # Bond information storage for click detection
        self._bonds: list[tuple[int, int, float]] = []  # (atom1_idx, atom2_idx, length)
        self._highlighted_bond: tuple[int, int] | None = None  # Currently clicked bond

        # Measurement mode
        self._measurement_mode = False
        self._measurement_points: list[tuple[float, float, float]] = []

        # Connect mouse events
        self.mpl_connect('button_press_event', self._on_click)
        self.mpl_connect('scroll_event', self._on_scroll)
        self.mpl_connect('motion_notify_event', self._on_motion)

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

        # Equal aspect ratio to prevent stretching
        self.ax.set_box_aspect([1, 1, 1])

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
        """Draw all slices in 3D.

        Always uses Z-slices through the density cube. The molecule rotation
        (for custom plane selection) is handled upstream when the density cube
        is calculated, not here.
        """
        self._setup_3d_axes()

        if not self._density_cube:
            self.ax.text2D(0.5, 0.5, 'No density data',
                          ha='center', va='center', transform=self.ax.transAxes,
                          fontsize=12, color='#888')
            self.draw()
            return

        # Always use Z-slices - molecule rotation is handled in density calculation
        slices = self._density_cube.get_z_slices(self._n_slices)

        # Handle empty slices
        if not slices:
            self.ax.text2D(0.5, 0.5, 'No slices available',
                          ha='center', va='center', transform=self.ax.transAxes,
                          fontsize=12, color='#888')
            self.draw()
            return

        # Get grid extent in Angstrom
        origin = self._density_cube.origin * BOHR_TO_ANGSTROM
        nx, ny, nz = self._density_cube.shape
        axes = self._density_cube.axes * BOHR_TO_ANGSTROM

        x_extent = nx * axes[0, 0]
        y_extent = ny * axes[1, 1]

        # Determine contour levels based on density type
        # Scale factor: <1 = tighter (more contours, lower values), >1 = looser (fewer, higher values)
        scale = self._contour_scale
        all_data = np.concatenate([s.flatten() for _, s in slices])

        if self._density_type == "deformation":
            max_abs = max(abs(all_data.min()), abs(all_data.max()))
            if max_abs > 0.001:
                # Base levels scaled - scale affects the threshold values
                base_levels = np.array([0.005, 0.01, 0.02, 0.04, 0.08])
                levels_pos = list(base_levels * scale)
                levels_neg = list(-base_levels[::-1] * scale)
            else:
                levels_pos = []
                levels_neg = []
        else:
            # Promolecule density
            max_val = all_data.max()
            if max_val > 0.001:
                # Generate contours from min to max, scaled
                # More contours when scale < 1, fewer when scale > 1
                n_contours = max(3, int(7 / scale))
                # Start from a fraction of max, go up to 80% of max
                start = max_val * 0.02 * scale
                end = max_val * 0.8
                if start < end:
                    levels_pos = list(np.linspace(start, end, n_contours))
                else:
                    levels_pos = [max_val * 0.5]
            else:
                levels_pos = []
            levels_neg = []

        # Draw slices
        self._draw_z_slices(slices, origin, x_extent, y_extent, levels_pos, levels_neg)

        # Calculate display limits
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

        # Calculate global min/max for consistent heatmap normalization
        if self._render_mode == "heatmap":
            all_data = np.concatenate([s.flatten() for _, s in slices])
            if self._density_type == "deformation":
                # For deformation, use symmetric normalization around zero
                max_abs = max(abs(all_data.min()), abs(all_data.max()))
                vmin, vmax = -max_abs, max_abs
            else:
                vmin, vmax = all_data.min(), all_data.max()
            norm = Normalize(vmin=vmin, vmax=vmax)
            cmap = cm.get_cmap(self._heatmap_colormap)

        for i, (z_coord, slice_data) in enumerate(slices):
            is_highlighted = (i == self._highlighted_slice)

            # Skip non-highlighted slices if "show only highlighted" is enabled
            if self._show_only_highlighted and not is_highlighted:
                continue

            alpha = 0.9 if is_highlighted else 0.15
            Z = np.full_like(X, z_coord)

            if self._render_mode == "heatmap":
                self._draw_slice_heatmap(X, Y, z_coord, slice_data, alpha, is_highlighted, norm, cmap)
            elif self._show_contours and len(levels_pos) > 0:
                self._draw_slice_contours(
                    X, Y, Z, slice_data, z_coord,
                    levels_pos, levels_neg,
                    alpha, is_highlighted
                )
            else:
                self._draw_slice_plane(X, Y, z_coord, slice_data, alpha, is_highlighted)

    def _draw_slice_contours(self, X, Y, Z, slice_data, z_coord,
                             levels_pos, levels_neg, alpha, is_highlighted):
        """Draw contour lines on a slice at given z-coordinate."""
        # First draw the plane behind the contours
        self._draw_slice_plane(X, Y, z_coord, slice_data, alpha * 0.3, is_highlighted)

        # For 3D contour with zdir='z', we pass X, Y grids and data directly
        # The data shape should match (nx, ny) where X is (nx, ny) and Y is (nx, ny)

        # Draw positive contours
        if levels_pos and len(levels_pos) > 0:
            valid_levels = [l for l in levels_pos if l < slice_data.max()]
            if valid_levels:
                try:
                    self.ax.contour(
                        X, Y, slice_data,
                        levels=valid_levels,
                        zdir='z',
                        offset=z_coord,
                        colors='#000000' if self._color_mode == 'bw' else '#2563eb',
                        linewidths=2.0 if is_highlighted else 0.6,
                        alpha=1.0 if is_highlighted else 0.4
                    )
                except Exception as e:
                    print(f"Contour error (pos): {e}")

        # Draw negative contours (for deformation density)
        if levels_neg and len(levels_neg) > 0:
            valid_levels = [l for l in levels_neg if l > slice_data.min()]
            if valid_levels:
                try:
                    self.ax.contour(
                        X, Y, slice_data,
                        levels=valid_levels,
                        zdir='z',
                        offset=z_coord,
                        colors='#666666' if self._color_mode == 'bw' else '#dc2626',
                        linestyles='dashed',
                        linewidths=2.0 if is_highlighted else 0.6,
                        alpha=1.0 if is_highlighted else 0.4
                    )
                except Exception as e:
                    print(f"Contour error (neg): {e}")

    def _draw_slice_plane(self, X, Y, z_coord, slice_data, alpha, is_highlighted):
        """Draw a semi-transparent slice plane as a simple rectangle."""
        # Get bounds
        x_min, x_max = X.min(), X.max()
        y_min, y_max = Y.min(), Y.max()

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

    def _draw_slice_heatmap(self, X, Y, z_coord, slice_data, alpha, is_highlighted,
                            norm, cmap):
        """Draw a slice as a color-mapped heatmap using batched quad patches.

        For 3D rendering, we divide the slice into small quads, each colored
        according to its density value. All quads are batched into a single
        Poly3DCollection for efficiency.
        """
        nx, ny = slice_data.shape

        # Subsample for performance - use every nth point
        step = max(1, min(nx, ny) // 20)  # Aim for ~20x20 quads max

        # Collect all quads and colors
        all_verts = []
        all_colors = []

        for i in range(0, nx - step, step):
            for j in range(0, ny - step, step):
                # Get corner coordinates
                x0, x1 = X[i, j], X[min(i + step, nx - 1), j]
                y0, y1 = Y[i, j], Y[i, min(j + step, ny - 1)]

                # Average density in this quad
                avg_density = slice_data[i:i+step, j:j+step].mean()

                # Get color from colormap
                color = cmap(norm(avg_density))

                # Add quad vertices
                all_verts.append([
                    (x0, y0, z_coord),
                    (x1, y0, z_coord),
                    (x1, y1, z_coord),
                    (x0, y1, z_coord)
                ])
                all_colors.append(color[:3])  # RGB only

        # Create single batched collection if we have quads
        if all_verts:
            face_alpha = 0.85 if is_highlighted else 0.25
            poly = Poly3DCollection(
                all_verts,
                alpha=face_alpha,
                facecolors=all_colors,
                edgecolors='none',
                linewidths=0
            )
            self.ax.add_collection3d(poly)

        # Draw edge around the highlighted slice using a line plot
        if is_highlighted:
            x_min, x_max = X.min(), X.max()
            y_min, y_max = Y.min(), Y.max()
            # Draw rectangle as lines
            self.ax.plot(
                [x_min, x_max, x_max, x_min, x_min],
                [y_min, y_min, y_max, y_max, y_min],
                [z_coord, z_coord, z_coord, z_coord, z_coord],
                color='#333333',
                linewidth=1.5,
                zorder=10
            )

    def _draw_bonds(self):
        """Draw bonds between atoms in 3D and store bond info for click detection."""
        self._bonds = []  # Clear previous bonds

        if not self._density_cube or not self._density_cube.atoms:
            return

        if not self._show_bonds:
            return

        atoms = self._density_cube.atoms
        n_atoms = len(atoms)

        for i in range(n_atoms):
            z1, x1, y1, z1_coord = atoms[i]
            x1 *= BOHR_TO_ANGSTROM
            y1 *= BOHR_TO_ANGSTROM
            z1_coord *= BOHR_TO_ANGSTROM
            r1 = COVALENT_RADII.get(z1, 1.5)

            for j in range(i + 1, n_atoms):
                z2, x2, y2, z2_coord = atoms[j]
                x2 *= BOHR_TO_ANGSTROM
                y2 *= BOHR_TO_ANGSTROM
                z2_coord *= BOHR_TO_ANGSTROM
                r2 = COVALENT_RADII.get(z2, 1.5)

                dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2_coord - z1_coord)**2)
                max_bond_dist = (r1 + r2) * 1.3

                if dist < max_bond_dist:
                    # Store bond info
                    self._bonds.append((i, j, dist))

                    # Check if this bond is highlighted
                    is_highlighted = (
                        self._highlighted_bond is not None and
                        ((i, j) == self._highlighted_bond or (j, i) == self._highlighted_bond)
                    )

                    if is_highlighted:
                        # Highlighted bond: thicker, colored
                        self.ax.plot(
                            [x1, x2], [y1, y2], [z1_coord, z2_coord],
                            color='#FF6600',  # Orange for highlighted
                            linewidth=4,
                            alpha=0.9,
                            zorder=5
                        )
                        # Show bond midpoint with density annotation
                        mid_x, mid_y, mid_z = (x1+x2)/2, (y1+y2)/2, (z1_coord+z2_coord)/2
                        if self._density_cube:
                            density = self._density_cube.get_density_at_point(mid_x, mid_y, mid_z)
                            label = f"{dist:.2f}Å\nρ={density:.3f}"
                            self.ax.text(mid_x, mid_y, mid_z + 0.3, label,
                                        fontsize=8, ha='center', color='#FF6600',
                                        fontweight='bold', zorder=20)
                    else:
                        # Normal bond
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
            x = ax * BOHR_TO_ANGSTROM
            y = ay * BOHR_TO_ANGSTROM
            z = az * BOHR_TO_ANGSTROM

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
        """Handle mouse click for atom/bond picking in 3D.

        Note: We use 'button_press_event' but only handle atom/bond picking.
        Matplotlib's default 3D rotation (drag) still works because we
        don't block the event - we just emit a signal if something is clicked.
        """
        # Only process left button clicks for picking
        if event.button != 1:  # 1 = left button
            return

        if event.inaxes != self.ax:
            return

        # Get click position in display coordinates
        click_x, click_y = event.x, event.y
        if click_x is None or click_y is None:
            return

        from mpl_toolkits.mplot3d import proj3d

        # Check for atom click first (if pick mode enabled)
        if self._pick_enabled and self._atom_positions_3d:
            min_dist = float('inf')
            closest_idx = -1

            for idx, (ax, ay, az) in enumerate(self._atom_positions_3d):
                x2d, y2d, _ = proj3d.proj_transform(ax, ay, az, self.ax.get_proj())
                x_disp, y_disp = self.ax.transData.transform((x2d, y2d))
                dist = np.sqrt((x_disp - click_x) ** 2 + (y_disp - click_y) ** 2)
                if dist < min_dist:
                    min_dist = dist
                    closest_idx = idx

            if closest_idx >= 0 and min_dist < 30:
                self.atom_picked.emit(closest_idx)
                return

        # Check for bond click (always enabled when bonds are shown)
        if self._show_bonds and self._bonds and self._atom_positions_3d:
            min_dist = float('inf')
            closest_bond = None

            for i, j, bond_length in self._bonds:
                if i >= len(self._atom_positions_3d) or j >= len(self._atom_positions_3d):
                    continue

                # Get bond midpoint
                ax1, ay1, az1 = self._atom_positions_3d[i]
                ax2, ay2, az2 = self._atom_positions_3d[j]
                mid_x, mid_y, mid_z = (ax1+ax2)/2, (ay1+ay2)/2, (az1+az2)/2

                # Project to display coordinates
                x2d, y2d, _ = proj3d.proj_transform(mid_x, mid_y, mid_z, self.ax.get_proj())
                x_disp, y_disp = self.ax.transData.transform((x2d, y2d))
                dist = np.sqrt((x_disp - click_x) ** 2 + (y_disp - click_y) ** 2)

                if dist < min_dist:
                    min_dist = dist
                    closest_bond = (i, j, bond_length)

            if closest_bond and min_dist < 40:  # 40 pixels threshold for bonds
                i, j, bond_length = closest_bond
                # Toggle highlight
                if self._highlighted_bond == (i, j) or self._highlighted_bond == (j, i):
                    self._highlighted_bond = None
                else:
                    self._highlighted_bond = (i, j)
                self._draw_slices()
                self.bond_clicked.emit(i, j, bond_length)

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

    def _on_motion(self, event):
        """Handle mouse motion for density hover display."""
        if not self._show_density_on_hover:
            return

        if event.inaxes != self.ax or not self._density_cube:
            if self._density_text:
                self._density_text.set_visible(False)
                self.draw_idle()
            return

        # For 3D, we can't easily get density under cursor
        # Instead, show density info in a fixed position
        # This requires tracking the current slice data
        # For now, we'll display general stats
        if self._density_text is None:
            self._density_text = self.fig.text(
                0.02, 0.98, '',
                fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
            )

        # Show density range info
        # Density is in e/Bohr³, convert to e/Å³ for display
        data = self._density_cube.data
        # Convert density from e/Bohr³ to e/Å³ (multiply by Bohr³/Å³)
        density_ang = data * (1.0 / BOHR_TO_ANGSTROM**3)
        # Calculate volume element for total electron count
        dV = abs(np.linalg.det(self._density_cube.axes))  # Volume in Bohr³
        total_electrons = data.sum() * dV
        info = (
            f"Density range:\n"
            f"  min: {density_ang.min():.4f} e/Å³\n"
            f"  max: {density_ang.max():.4f} e/Å³\n"
            f"  total: {total_electrons:.1f} e"
        )
        self._density_text.set_text(info)
        self._density_text.set_visible(True)
        self.draw_idle()

    def set_show_density_hover(self, show: bool):
        """Enable or disable density display on hover."""
        self._show_density_on_hover = show
        if not show and self._density_text:
            self._density_text.set_visible(False)
            self.draw_idle()

    def set_contour_scale(self, scale: float):
        """Set contour level scale (0.25-4.0).

        Args:
            scale: <1 = tighter/more contours, >1 = looser/fewer contours
        """
        self._contour_scale = max(0.25, min(4.0, scale))
        if self._density_cube:
            self._draw_slices()

    def set_render_mode(self, mode: str):
        """Set render mode ('contour' or 'heatmap').

        Args:
            mode: 'contour' for contour lines, 'heatmap' for color-mapped bitmap
        """
        self._render_mode = mode
        if self._density_cube:
            self._draw_slices()

    def set_heatmap_colormap(self, colormap: str):
        """Set the colormap for heatmap mode.

        Args:
            colormap: matplotlib colormap name (e.g., 'RdYlBu_r', 'viridis', 'hot')
        """
        self._heatmap_colormap = colormap
        if self._density_cube and self._render_mode == "heatmap":
            self._draw_slices()

    def get_render_mode(self) -> str:
        """Get current render mode."""
        return self._render_mode

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

        # Maintain equal aspect ratio based on data ranges
        x_range = 2 * half_x
        y_range = 2 * half_y
        z_range = 2 * half_z
        max_range = max(x_range, y_range, z_range)
        if max_range > 0:
            self.ax.set_box_aspect([x_range / max_range,
                                    y_range / max_range,
                                    z_range / max_range])

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

    def reset_view(self):
        """Reset zoom and rotation to defaults."""
        self._zoom_level = 1.0
        self.ax.view_init(elev=25, azim=-60)
        self._apply_zoom()

    def get_slice_z_coord(self, slice_index: int) -> float | None:
        """Get the Z-coordinate (in Å) for a given slice index."""
        if not self._density_cube:
            return None
        slices = self._density_cube.get_z_slices(self._n_slices)
        if not slices or slice_index >= len(slices):
            return None
        return slices[slice_index][0]

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

    def set_custom_plane(self, plane_definition):
        """Legacy method - custom plane is now handled by recalculating the density cube.

        This method is kept for API compatibility but does nothing.
        The plane rotation is applied when the density cube is calculated,
        so slicing always uses simple Z-slices.
        """
        # No-op: plane rotation is handled upstream in density calculation
        pass

    def set_show_bonds(self, show: bool):
        """Enable or disable bond visualization."""
        self._show_bonds = show
        if self._density_cube:
            self._draw_slices()

    def set_show_only_highlighted(self, show_only: bool):
        """Toggle showing only the highlighted slice vs all slices."""
        self._show_only_highlighted = show_only
        if self._density_cube:
            self._draw_slices()

    def get_current_slice_data(self) -> tuple[np.ndarray, float] | None:
        """Get the data for the currently highlighted slice.

        Returns:
            (slice_data, z_coord) or None if no data
        """
        if not self._density_cube:
            return None
        slices = self._density_cube.get_z_slices(self._n_slices)
        if not slices or self._highlighted_slice >= len(slices):
            return None
        return slices[self._highlighted_slice]

    def get_density_type(self) -> str:
        """Get current density type."""
        return self._density_type

    def get_color_mode(self) -> str:
        """Get current color mode."""
        return self._color_mode

    def get_contour_scale(self) -> float:
        """Get current contour scale."""
        return self._contour_scale

    def get_density_cube(self) -> DensityCube | None:
        """Get the current density cube."""
        return self._density_cube

    def clear_bond_highlight(self):
        """Clear any highlighted bond."""
        if self._highlighted_bond is not None:
            self._highlighted_bond = None
            if self._density_cube:
                self._draw_slices()

    def get_bond_info(self, atom1_idx: int, atom2_idx: int) -> dict | None:
        """Get detailed info about a bond.

        Returns dict with: length, midpoint, density_at_midpoint, atom1_symbol, atom2_symbol
        """
        if not self._density_cube or not self._density_cube.atoms:
            return None

        atoms = self._density_cube.atoms
        if atom1_idx >= len(atoms) or atom2_idx >= len(atoms):
            return None

        z1, x1, y1, z1c = atoms[atom1_idx]
        z2, x2, y2, z2c = atoms[atom2_idx]

        x1, y1, z1c = x1 * BOHR_TO_ANGSTROM, y1 * BOHR_TO_ANGSTROM, z1c * BOHR_TO_ANGSTROM
        x2, y2, z2c = x2 * BOHR_TO_ANGSTROM, y2 * BOHR_TO_ANGSTROM, z2c * BOHR_TO_ANGSTROM

        length = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2c-z1c)**2)
        mid_x, mid_y, mid_z = (x1+x2)/2, (y1+y2)/2, (z1c+z2c)/2
        density = self._density_cube.get_density_at_point(mid_x, mid_y, mid_z)

        return {
            'length': length,
            'midpoint': (mid_x, mid_y, mid_z),
            'density': density,
            'atom1_symbol': Z_TO_SYMBOL.get(z1, '?'),
            'atom2_symbol': Z_TO_SYMBOL.get(z2, '?'),
            'atom1_z': z1,
            'atom2_z': z2,
        }


class SliceExplorer(QWidget):
    """Interactive 3D Slice Explorer widget."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._promolecule_cube: DensityCube | None = None
        self._deformation_cube: DensityCube | None = None

        # Animation state
        self._animation_timer = QTimer(self)
        self._animation_timer.timeout.connect(self._animation_step)
        self._animation_progress = 0.0
        self._animation_direction = 1  # 1 = forward (pro->def), -1 = backward
        self._is_animating = False

        # 2D preview dialog (lazy init)
        self._preview_dialog: Slice2DPreviewDialog | None = None

        self._setup_ui()

    def _setup_ui(self):
        """Set up the widget UI with a cleaner two-row layout."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        # 3D Canvas
        self.canvas = SliceExplorerCanvas(self)
        self.canvas.setMinimumSize(400, 400)
        layout.addWidget(self.canvas, stretch=1)

        # === ROW 1: View Options ===
        row1 = QHBoxLayout()
        row1.setSpacing(20)

        # Density type selector
        type_group = QHBoxLayout()
        type_group.setSpacing(6)
        type_label = QLabel("Density:")
        type_label.setStyleSheet("font-weight: bold;")
        type_label.setToolTip("Switch between promolecule and deformation density")
        type_group.addWidget(type_label)
        self.density_combo = QComboBox()
        self.density_combo.addItems(["Promolecule", "Deformation"])
        self.density_combo.setToolTip("Promolecule: atom positions | Deformation: bonding effects")
        self.density_combo.currentIndexChanged.connect(self._on_density_type_changed)
        type_group.addWidget(self.density_combo)
        row1.addLayout(type_group)

        # Render mode toggle buttons (Contour vs Heatmap)
        render_group = QHBoxLayout()
        render_group.setSpacing(2)
        self.contour_mode_btn = QPushButton("Contour")
        self.contour_mode_btn.setCheckable(True)
        self.contour_mode_btn.setChecked(True)
        self.contour_mode_btn.setToolTip("Show contour lines (good for printing)")
        self.contour_mode_btn.setMaximumWidth(60)
        self.heatmap_mode_btn = QPushButton("Heatmap")
        self.heatmap_mode_btn.setCheckable(True)
        self.heatmap_mode_btn.setToolTip("Show color gradient (good for screen visualization)")
        self.heatmap_mode_btn.setMaximumWidth(60)
        # Button group for exclusive selection
        self.render_mode_group = QButtonGroup(self)
        self.render_mode_group.addButton(self.contour_mode_btn, 0)
        self.render_mode_group.addButton(self.heatmap_mode_btn, 1)
        self.render_mode_group.idClicked.connect(self._on_render_mode_changed)
        render_group.addWidget(self.contour_mode_btn)
        render_group.addWidget(self.heatmap_mode_btn)
        row1.addLayout(render_group)

        # Color mode toggle buttons (only for contour mode)
        color_group = QHBoxLayout()
        color_group.setSpacing(2)
        self.bw_btn = QPushButton("B&&W")
        self.bw_btn.setCheckable(True)
        self.bw_btn.setChecked(True)
        self.bw_btn.setToolTip("Black and white contours")
        self.bw_btn.setMaximumWidth(45)
        self.color_btn = QPushButton("Color")
        self.color_btn.setCheckable(True)
        self.color_btn.setToolTip("Colored contours (blue=positive, red=negative)")
        self.color_btn.setMaximumWidth(50)
        # Button group for exclusive selection
        self.color_button_group = QButtonGroup(self)
        self.color_button_group.addButton(self.bw_btn, 0)
        self.color_button_group.addButton(self.color_btn, 1)
        self.color_button_group.idClicked.connect(self._on_color_mode_changed)
        color_group.addWidget(self.bw_btn)
        color_group.addWidget(self.color_btn)
        row1.addLayout(color_group)

        # Separator
        sep1 = QFrame()
        sep1.setFrameShape(QFrame.Shape.VLine)
        sep1.setFrameShadow(QFrame.Shadow.Sunken)
        row1.addWidget(sep1)

        # Display options (checkboxes)
        self.contours_check = QCheckBox("Contours")
        self.contours_check.setChecked(True)
        self.contours_check.setToolTip("Show/hide contour lines")
        self.contours_check.stateChanged.connect(self._on_contours_toggled)
        row1.addWidget(self.contours_check)

        self.bonds_check = QCheckBox("Bonds")
        self.bonds_check.setChecked(True)
        self.bonds_check.setToolTip("Show/hide molecular bonds")
        self.bonds_check.stateChanged.connect(self._on_bonds_toggled)
        row1.addWidget(self.bonds_check)

        self.single_slice_check = QCheckBox("Single")
        self.single_slice_check.setChecked(False)
        self.single_slice_check.setToolTip("Show only the selected slice, hide all others")
        self.single_slice_check.stateChanged.connect(self._on_single_slice_toggled)
        row1.addWidget(self.single_slice_check)

        self.density_info_check = QCheckBox("Info")
        self.density_info_check.setChecked(False)
        self.density_info_check.setToolTip("Show electron density statistics")
        self.density_info_check.stateChanged.connect(self._on_density_info_toggled)
        row1.addWidget(self.density_info_check)

        # Separator
        sep2 = QFrame()
        sep2.setFrameShape(QFrame.Shape.VLine)
        sep2.setFrameShadow(QFrame.Shadow.Sunken)
        row1.addWidget(sep2)

        # Contour spacing slider
        spacing_group = QHBoxLayout()
        spacing_group.setSpacing(4)
        spacing_label = QLabel("Spacing:")
        spacing_label.setToolTip("Adjust contour level spacing")
        spacing_group.addWidget(spacing_label)
        self.contour_scale_slider = QSlider(Qt.Orientation.Horizontal)
        self.contour_scale_slider.setMinimum(25)  # 0.25x
        self.contour_scale_slider.setMaximum(200)  # 2.0x
        self.contour_scale_slider.setValue(100)  # 1.0x (default)
        self.contour_scale_slider.setFixedWidth(80)
        self.contour_scale_slider.setToolTip("Left = more contours (tighter), Right = fewer contours (looser)")
        self.contour_scale_slider.valueChanged.connect(self._on_contour_scale_changed)
        spacing_group.addWidget(self.contour_scale_slider)
        self.contour_scale_label = QLabel("1.0x")
        self.contour_scale_label.setMinimumWidth(30)
        spacing_group.addWidget(self.contour_scale_label)
        row1.addLayout(spacing_group)

        row1.addStretch()
        layout.addLayout(row1)

        # === ROW 2: Navigation Controls ===
        row2 = QHBoxLayout()
        row2.setSpacing(20)

        # Slice navigation
        slice_group = QHBoxLayout()
        slice_group.setSpacing(6)
        slice_label = QLabel("Slice:")
        slice_label.setStyleSheet("font-weight: bold;")
        slice_label.setToolTip("Navigate through density slices")
        slice_group.addWidget(slice_label)

        self.slice_slider = QSlider(Qt.Orientation.Horizontal)
        self.slice_slider.setMinimum(0)
        self.slice_slider.setMaximum(14)
        self.slice_slider.setValue(7)
        self.slice_slider.setMinimumWidth(150)
        self.slice_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.slice_slider.setTickInterval(1)
        self.slice_slider.setToolTip("Drag to select slice")
        self.slice_slider.valueChanged.connect(self._on_slice_changed)
        slice_group.addWidget(self.slice_slider, stretch=1)

        self.slice_label = QLabel("8/15")
        self.slice_label.setMinimumWidth(90)
        self.slice_label.setToolTip("Current slice / total slices (Z-coordinate)")
        slice_group.addWidget(self.slice_label)
        row2.addLayout(slice_group, stretch=1)

        # Separator
        sep3 = QFrame()
        sep3.setFrameShape(QFrame.Shape.VLine)
        sep3.setFrameShadow(QFrame.Shadow.Sunken)
        row2.addWidget(sep3)

        # Zoom controls
        zoom_group = QHBoxLayout()
        zoom_group.setSpacing(2)
        zoom_label = QLabel("Zoom:")
        zoom_label.setToolTip("Zoom in/out (or use scroll wheel)")
        zoom_group.addWidget(zoom_label)

        zoom_in_btn = QPushButton("+")
        zoom_in_btn.setFixedWidth(28)
        zoom_in_btn.setToolTip("Zoom in")
        zoom_in_btn.clicked.connect(self.canvas.zoom_in)
        zoom_group.addWidget(zoom_in_btn)

        zoom_out_btn = QPushButton("-")
        zoom_out_btn.setFixedWidth(28)
        zoom_out_btn.setToolTip("Zoom out")
        zoom_out_btn.clicked.connect(self.canvas.zoom_out)
        zoom_group.addWidget(zoom_out_btn)
        row2.addLayout(zoom_group)

        # Rotation controls
        rotate_group = QHBoxLayout()
        rotate_group.setSpacing(2)
        rotate_label = QLabel("Rotate:")
        rotate_label.setToolTip("Rotate view (or drag mouse on canvas)")
        rotate_group.addWidget(rotate_label)

        rotate_left_btn = QPushButton("\u2190")
        rotate_left_btn.setFixedWidth(28)
        rotate_left_btn.setToolTip("Rotate left")
        rotate_left_btn.clicked.connect(self.canvas.rotate_left)
        rotate_group.addWidget(rotate_left_btn)

        rotate_right_btn = QPushButton("\u2192")
        rotate_right_btn.setFixedWidth(28)
        rotate_right_btn.setToolTip("Rotate right")
        rotate_right_btn.clicked.connect(self.canvas.rotate_right)
        rotate_group.addWidget(rotate_right_btn)

        rotate_up_btn = QPushButton("\u2191")
        rotate_up_btn.setFixedWidth(28)
        rotate_up_btn.setToolTip("Rotate up")
        rotate_up_btn.clicked.connect(self.canvas.rotate_up)
        rotate_group.addWidget(rotate_up_btn)

        rotate_down_btn = QPushButton("\u2193")
        rotate_down_btn.setFixedWidth(28)
        rotate_down_btn.setToolTip("Rotate down")
        rotate_down_btn.clicked.connect(self.canvas.rotate_down)
        rotate_group.addWidget(rotate_down_btn)
        row2.addLayout(rotate_group)

        # Reset View button
        reset_btn = QPushButton("Reset View")
        reset_btn.setToolTip("Reset zoom and rotation to default")
        reset_btn.clicked.connect(self.canvas.reset_view)
        row2.addWidget(reset_btn)

        # Separator
        sep4 = QFrame()
        sep4.setFrameShape(QFrame.Shape.VLine)
        sep4.setFrameShadow(QFrame.Shadow.Sunken)
        row2.addWidget(sep4)

        # View 2D button
        self.view_2d_btn = QPushButton("View 2D")
        self.view_2d_btn.setToolTip("Show current slice as a 2D contour plot")
        self.view_2d_btn.clicked.connect(self._show_2d_preview)
        row2.addWidget(self.view_2d_btn)

        # Animate button (morph between promolecule and deformation)
        self.animate_btn = QPushButton("Animate")
        self.animate_btn.setToolTip("Animate transition between promolecule and deformation density")
        self.animate_btn.setCheckable(True)
        self.animate_btn.clicked.connect(self._toggle_animation)
        row2.addWidget(self.animate_btn)

        layout.addLayout(row2)

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

    def _on_color_mode_changed(self, button_id: int):
        """Handle color mode toggle button change."""
        mode = "bw" if button_id == 0 else "color"
        self.canvas.set_color_mode(mode)

    def _on_render_mode_changed(self, button_id: int):
        """Handle render mode toggle button change."""
        mode = "contour" if button_id == 0 else "heatmap"
        self.canvas.set_render_mode(mode)
        # Disable/enable contour-specific controls
        is_contour = (mode == "contour")
        self.bw_btn.setEnabled(is_contour)
        self.color_btn.setEnabled(is_contour)
        self.contours_check.setEnabled(is_contour)
        self.contour_scale_slider.setEnabled(is_contour)

    def _on_contours_toggled(self, state: int):
        """Handle contours checkbox toggle."""
        self.canvas.set_show_contours(state == Qt.CheckState.Checked.value)

    def _on_density_info_toggled(self, state: int):
        """Handle density info checkbox toggle."""
        self.canvas.set_show_density_hover(state == Qt.CheckState.Checked.value)

    def _on_contour_scale_changed(self, value: int):
        """Handle contour scale slider change."""
        scale = value / 100.0
        self.contour_scale_label.setText(f"{scale:.1f}x")
        self.canvas.set_contour_scale(scale)

    def _on_bonds_toggled(self, state: int):
        """Handle bonds checkbox toggle."""
        self.canvas.set_show_bonds(state == Qt.CheckState.Checked.value)

    def _on_single_slice_toggled(self, state: int):
        """Handle single slice checkbox toggle."""
        self.canvas.set_show_only_highlighted(state == Qt.CheckState.Checked.value)

    def _update_slice_label(self):
        """Update the slice number label with Z-coordinate."""
        current = self.slice_slider.value() + 1
        total = self.slice_slider.maximum() + 1
        z_coord = self.canvas.get_slice_z_coord(self.slice_slider.value())
        if z_coord is not None:
            self.slice_label.setText(f"{current}/{total} (z={z_coord:.1f}Å)")
        else:
            self.slice_label.setText(f"{current}/{total}")

    def clear(self):
        """Clear the explorer."""
        self._promolecule_cube = None
        self._deformation_cube = None
        self.canvas._density_cube = None
        self.canvas._setup_3d_axes()
        self.canvas.draw()
        self._stop_animation()

    def _show_2d_preview(self):
        """Show the current slice in a 2D preview dialog."""
        slice_result = self.canvas.get_current_slice_data()
        if slice_result is None:
            return

        z_coord, slice_data = slice_result
        cube = self.canvas.get_density_cube()

        if self._preview_dialog is None:
            self._preview_dialog = Slice2DPreviewDialog(self)

        self._preview_dialog.show_slice(
            slice_data=slice_data,
            z_coord=z_coord,
            density_type=self.canvas.get_density_type(),
            color_mode=self.canvas.get_color_mode(),
            atoms=cube.atoms if cube else None,
            contour_scale=self.canvas.get_contour_scale(),
            origin=cube.origin if cube else None,
            axes=cube.axes if cube else None,
            render_mode=self.canvas.get_render_mode(),
            colormap=self.canvas._heatmap_colormap,
        )
        self._preview_dialog.show()
        self._preview_dialog.raise_()

    def _toggle_animation(self, checked: bool):
        """Start or stop the density animation."""
        if checked:
            self._start_animation()
        else:
            self._stop_animation()

    def _start_animation(self):
        """Start animating between promolecule and deformation density."""
        if not self._promolecule_cube or not self._deformation_cube:
            self.animate_btn.setChecked(False)
            return

        self._is_animating = True
        self._animation_progress = 0.0
        self._animation_direction = 1
        self._animation_timer.start(100)  # 100ms per frame = 10 fps

    def _stop_animation(self):
        """Stop the density animation."""
        self._animation_timer.stop()
        self._is_animating = False
        self.animate_btn.setChecked(False)

    def _animation_step(self):
        """One step of the animation - blend between densities."""
        if not self._promolecule_cube or not self._deformation_cube:
            self._stop_animation()
            return

        # Update progress
        self._animation_progress += 0.05 * self._animation_direction

        # Reverse at boundaries
        if self._animation_progress >= 1.0:
            self._animation_progress = 1.0
            self._animation_direction = -1
        elif self._animation_progress <= 0.0:
            self._animation_progress = 0.0
            self._animation_direction = 1

        # Create blended density cube
        t = self._animation_progress
        # Blend: at t=0 show promolecule, at t=1 show molecular (promolecule + deformation)
        # This shows the "formation of bonds" effect
        blended_data = self._promolecule_cube.data + t * self._deformation_cube.data

        blended_cube = DensityCube(
            origin=self._promolecule_cube.origin,
            axes=self._promolecule_cube.axes,
            data=blended_data,
            atoms=self._promolecule_cube.atoms,
        )

        # Update display
        density_type = "deformation" if t > 0.5 else "promolecule"
        self.canvas.set_density_cube(blended_cube, density_type)

        # Update combo to show progress
        self.density_combo.blockSignals(True)
        if t < 0.5:
            self.density_combo.setCurrentText("Promolecule")
        else:
            # Try to select Deformation if available
            idx = self.density_combo.findText("Deformation")
            if idx >= 0:
                self.density_combo.setCurrentIndex(idx)
        self.density_combo.blockSignals(False)
