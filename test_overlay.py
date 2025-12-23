"""
Test script to generate overlaid promolecule + deformation contours.
Helps verify coordinate alignment between the two density calculations.
"""

import tempfile
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from dorothy.core.cif_parser import parse_cif, load_example_molecules
from dorothy.core.xtb_manager import is_xtb_installed, run_xtb_density
from dorothy.core.density import (
    parse_cube_file,
    calculate_deformation_density,
    create_density_cube_from_structure,
    calculate_promolecule_density,
)


def generate_overlay(structure, output_path: Path):
    """Generate PDF with overlaid promolecule (BW) and deformation (color) contours."""

    if not is_xtb_installed():
        print("xTB not installed - can't generate deformation density")
        return

    # Run xTB to get molecular density
    print("Running xTB calculation...")
    with tempfile.TemporaryDirectory() as tmp_dir:
        cube_path = run_xtb_density(structure, Path(tmp_dir))
        if not cube_path:
            print("xTB calculation failed")
            return

        molecular_density = parse_cube_file(cube_path)

    # Get coords from cube file
    z_to_symbol = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I'}
    symbols = [z_to_symbol.get(atom[0], 'C') for atom in molecular_density.atoms]
    coords = np.array([[atom[1], atom[2], atom[3]] for atom in molecular_density.atoms]) * 0.529177

    # Calculate promolecule on same grid
    print("Calculating promolecule on same grid...")
    promolecule_data = calculate_promolecule_density(
        coords,
        symbols,
        molecular_density.origin,
        molecular_density.axes,
        molecular_density.shape,
    )

    # Debug: Compare totals
    mol_sum = np.sum(molecular_density.data)
    pro_sum = np.sum(promolecule_data)
    print(f"\n=== Density comparison ===")
    print(f"Molecular density sum: {mol_sum:.4f}")
    print(f"Promolecule density sum: {pro_sum:.4f}")
    print(f"Ratio (mol/pro): {mol_sum/pro_sum:.4f}")
    print(f"Molecular max: {molecular_density.data.max():.6f}")
    print(f"Promolecule max: {promolecule_data.max():.6f}")
    print(f"Molecular min: {molecular_density.data.min():.6f}")
    print(f"Promolecule min: {promolecule_data.min():.6f}")

    # Calculate deformation WITHOUT normalization to see raw difference
    deform_raw = molecular_density.data - promolecule_data
    print(f"\nRaw deformation (no normalization):")
    print(f"  Max: {deform_raw.max():.6f}")
    print(f"  Min: {deform_raw.min():.6f}")
    print(f"  Sum: {deform_raw.sum():.6f}")

    # Normalized promolecule
    promolecule_norm = promolecule_data * (mol_sum / pro_sum)
    deform_norm = molecular_density.data - promolecule_norm
    print(f"\nNormalized deformation:")
    print(f"  Max: {deform_norm.max():.6f}")
    print(f"  Min: {deform_norm.min():.6f}")
    print(f"  Sum: {deform_norm.sum():.6f}")

    # Check where positive deformation occurs
    positive_mask = deform_norm > 0.01
    print(f"\nPositive deformation (>0.01): {positive_mask.sum()} grid points")
    print(f"Negative deformation (<-0.01): {(deform_norm < -0.01).sum()} grid points")

    # Sample a middle slice
    mid_z = molecular_density.shape[2] // 2
    mid_slice_mol = molecular_density.data[:, :, mid_z]
    mid_slice_pro = promolecule_norm[:, :, mid_z]
    mid_slice_def = deform_norm[:, :, mid_z]
    print(f"\nMiddle slice (z={mid_z}):")
    print(f"  Molecular max: {mid_slice_mol.max():.4f}")
    print(f"  Promolecule max: {mid_slice_pro.max():.4f}")
    print(f"  Deformation range: [{mid_slice_def.min():.4f}, {mid_slice_def.max():.4f}]")

    # Use NORMALIZED deformation for visualization
    deformation_data = deform_norm

    # Generate more slices for better resolution
    n_slices = 25  # 5x5 grid
    nz = molecular_density.shape[2]
    # Cover more of the z-range to catch π-density above/below plane
    slice_indices = np.linspace(nz // 6, 5 * nz // 6, n_slices, dtype=int)

    # Create PDF with multiple pages
    print(f"\nGenerating overlay PDF: {output_path}")
    # Fixed contour levels for consistency across all slices
    # Use logarithmic-like spacing to show both core and bonding regions
    # Cap the maximum to avoid nuclear peaks dominating
    DENSITY_CAP = 0.3  # Cap density display at this value
    DEFORM_CAP = 0.15  # Cap deformation display

    with PdfPages(output_path) as pdf:
        # Page 1: Deformation overlay (5x5 grid)
        fig, axes = plt.subplots(5, 5, figsize=(15, 15))
        fig.suptitle(f"{structure.name} - Deformation Density (blue=accumulation, red=depletion)", fontsize=14)

        # Fixed levels for deformation - same across all slices for comparability
        pos_levels = [0.01, 0.02, 0.04, 0.08, 0.12]
        neg_levels = [-0.12, -0.08, -0.04, -0.02, -0.01]

        for idx, (ax, slice_idx) in enumerate(zip(axes.flat, slice_indices)):
            deform_slice = deformation_data[:, :, slice_idx]

            z_bohr = molecular_density.origin[2] + slice_idx * molecular_density.axes[2, 2]
            z_ang = z_bohr * 0.529177

            # Deformation contours with fixed levels
            try:
                ax.contour(deform_slice.T, levels=pos_levels, colors='blue',
                          linewidths=0.8, linestyles='solid')
            except:
                pass
            try:
                ax.contour(deform_slice.T, levels=neg_levels, colors='red',
                          linewidths=0.8, linestyles='dashed')
            except:
                pass

            # Light promolecule outline for reference
            promol_slice = promolecule_norm[:, :, slice_idx]
            try:
                ax.contour(promol_slice.T, levels=[0.05, 0.1], colors='gray',
                          linewidths=0.3, linestyles='solid', alpha=0.5)
            except:
                pass

            # Atom positions
            for atom in molecular_density.atoms:
                z_atom = atom[3] * 0.529177
                if abs(z_atom - z_ang) < 0.3:
                    x_grid = (atom[1] - molecular_density.origin[0]) / molecular_density.axes[0, 0]
                    y_grid = (atom[2] - molecular_density.origin[1]) / molecular_density.axes[1, 1]
                    ax.plot(x_grid, y_grid, 'k.', markersize=2)

            ax.set_aspect('equal')
            ax.set_title(f"z={z_ang:.1f}Å", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()

        # Page 2: Raw molecular density with capped levels (5x5 grid)
        fig2, axes2 = plt.subplots(5, 5, figsize=(15, 15))
        fig2.suptitle(f"{structure.name} - xTB Molecular Density (capped at {DENSITY_CAP})", fontsize=14)

        # Fixed levels - don't go too high to avoid nuclear peaks
        mol_levels = [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25]

        for idx, (ax, slice_idx) in enumerate(zip(axes2.flat, slice_indices)):
            mol_slice = molecular_density.data[:, :, slice_idx]

            z_bohr = molecular_density.origin[2] + slice_idx * molecular_density.axes[2, 2]
            z_ang = z_bohr * 0.529177

            # Molecular density contours with capped levels
            try:
                ax.contour(mol_slice.T, levels=mol_levels, colors='black', linewidths=0.6)
            except:
                pass

            # Atom positions
            for atom in molecular_density.atoms:
                z_atom = atom[3] * 0.529177
                if abs(z_atom - z_ang) < 0.3:
                    x_grid = (atom[1] - molecular_density.origin[0]) / molecular_density.axes[0, 0]
                    y_grid = (atom[2] - molecular_density.origin[1]) / molecular_density.axes[1, 1]
                    ax.plot(x_grid, y_grid, 'r.', markersize=2)

            ax.set_aspect('equal')
            ax.set_title(f"z={z_ang:.1f}Å max={mol_slice.max():.2f}", fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])

        plt.tight_layout()
        pdf.savefig(fig2)
        plt.close()

    print(f"Done!")


if __name__ == "__main__":
    examples = load_example_molecules(Path(__file__).parent / "examples")

    if not examples:
        print("No example molecules found in examples/ directory")
        exit(1)

    structure = examples[0]
    print(f"Using molecule: {structure.name}")

    output_path = Path(__file__).parent / "test_output" / "overlay_test.pdf"
    output_path.parent.mkdir(exist_ok=True)

    generate_overlay(structure, output_path)
