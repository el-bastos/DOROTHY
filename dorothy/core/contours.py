"""
Contour slice generation and PDF export.

Creates electron density contour plots suitable for printing on transparent sheets.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

from dorothy.core.density import DensityCube
from dorothy.core.constants import (
    CONTOUR_LEVELS_PROMOLECULE,
    CONTOUR_LEVELS_DEFORMATION,
    DEFAULT_N_SLICES,
    PAGE_SIZE_A5,
)


@dataclass
class ContourSettings:
    """Settings for contour generation."""
    n_slices: int = DEFAULT_N_SLICES
    color_mode: str = "bw"  # "bw" or "color"
    n_contour_levels: int = 10
    page_size: tuple[float, float] = PAGE_SIZE_A5
    detail_level: str = "simple"  # "simple" or "advanced"
    # simple: fewer contours, auto-scaled per slice (easier for students)
    # advanced: fixed contour levels across all slices, shows π-bonds


def generate_promolecule_contours(
    density: DensityCube,
    output_dir: Path,
    settings: ContourSettings,
    progress_callback: Optional[callable] = None,
) -> list[Path]:
    """
    Generate contour plots for promolecule density.

    Args:
        density: Density cube
        output_dir: Output directory
        settings: Contour settings
        progress_callback: Optional callback(current, total, message)

    Returns:
        List of generated PDF paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    slices = density.get_z_slices(settings.n_slices)
    pdf_paths = []

    # Determine contour levels based on detail level
    if settings.detail_level == "advanced":
        # Fixed levels across all slices - capped to show bonding regions
        levels = CONTOUR_LEVELS_PROMOLECULE
    else:
        # Simple mode: auto-scaled based on data percentiles
        all_data = np.concatenate([s[1].flatten() for s in slices])
        positive_data = all_data[all_data > 0]
        if len(positive_data) > 0:
            vmin, vmax = np.percentile(positive_data, [5, 95])
            levels = np.linspace(vmin, vmax, settings.n_contour_levels)
        else:
            levels = np.linspace(0.01, 1.0, settings.n_contour_levels)

    for i, (z_coord, slice_data) in enumerate(slices):
        if progress_callback:
            progress_callback(i + 1, settings.n_slices, f"Generating slice {i + 1}/{settings.n_slices}")

        pdf_path = output_dir / f"slice_{i + 1:02d}.pdf"

        fig, ax = plt.subplots(figsize=settings.page_size)
        fig.patch.set_facecolor('white')

        # Create contour plot
        try:
            if settings.color_mode == "bw":
                ax.contour(slice_data.T, levels=levels, colors='black', linewidths=0.5)
            else:
                ax.contour(slice_data.T, levels=levels, cmap='Blues', linewidths=0.5)
        except ValueError:
            # No contours at these levels for this slice
            pass

        # Style
        ax.set_aspect('equal')
        ax.axis('off')

        # Add slice label
        ax.text(0.02, 0.98, f"Slice {i + 1}/{settings.n_slices}",
                transform=ax.transAxes, fontsize=8, va='top',
                color='gray')
        ax.text(0.02, 0.02, f"z = {z_coord:.2f} A",
                transform=ax.transAxes, fontsize=8, va='bottom',
                color='gray')

        plt.tight_layout(pad=0.5)
        fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
        plt.close(fig)

        pdf_paths.append(pdf_path)

    return pdf_paths


def generate_deformation_contours(
    density: DensityCube,
    output_dir: Path,
    settings: ContourSettings,
    progress_callback: Optional[callable] = None,
) -> list[Path]:
    """
    Generate contour plots for deformation density.

    Positive (accumulation) = solid lines / blue
    Negative (depletion) = dashed lines / red
    Zero contour omitted

    Args:
        density: Deformation density cube
        output_dir: Output directory
        settings: Contour settings
        progress_callback: Optional callback(current, total, message)

    Returns:
        List of generated PDF paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    slices = density.get_z_slices(settings.n_slices)
    pdf_paths = []

    # Determine contour levels based on detail level
    if settings.detail_level == "advanced":
        # Fixed levels - these show both σ and π bonding features
        positive_levels = CONTOUR_LEVELS_DEFORMATION
        negative_levels = [-l for l in reversed(CONTOUR_LEVELS_DEFORMATION)]
    else:
        # Simple mode: auto-scaled based on data, fewer levels
        all_data = np.concatenate([s[1].flatten() for s in slices])
        max_abs = np.percentile(np.abs(all_data), 95)
        n_levels = settings.n_contour_levels // 2
        positive_levels = np.linspace(max_abs * 0.1, max_abs, n_levels).tolist()
        negative_levels = (-np.array(positive_levels[::-1])).tolist()

    for i, (z_coord, slice_data) in enumerate(slices):
        if progress_callback:
            progress_callback(i + 1, settings.n_slices, f"Generating slice {i + 1}/{settings.n_slices}")

        pdf_path = output_dir / f"slice_{i + 1:02d}.pdf"

        fig, ax = plt.subplots(figsize=settings.page_size)
        fig.patch.set_facecolor('white')

        # Create contour plots
        try:
            if settings.color_mode == "bw":
                # Positive: solid black
                ax.contour(slice_data.T, levels=positive_levels,
                          colors='black', linewidths=0.5, linestyles='solid')
                # Negative: dashed black
                ax.contour(slice_data.T, levels=negative_levels,
                          colors='black', linewidths=0.5, linestyles='dashed')
            else:
                # Positive: blue
                ax.contour(slice_data.T, levels=positive_levels,
                          colors='blue', linewidths=0.5, linestyles='solid')
                # Negative: red
                ax.contour(slice_data.T, levels=negative_levels,
                          colors='red', linewidths=0.5, linestyles='dashed')
        except ValueError:
            # No contours at these levels for this slice
            pass

        # Style
        ax.set_aspect('equal')
        ax.axis('off')

        # Add slice label
        ax.text(0.02, 0.98, f"Slice {i + 1}/{settings.n_slices}",
                transform=ax.transAxes, fontsize=8, va='top',
                color='gray')
        ax.text(0.02, 0.02, f"z = {z_coord:.2f} A",
                transform=ax.transAxes, fontsize=8, va='bottom',
                color='gray')

        # Legend for deformation
        if settings.color_mode == "bw":
            legend_text = "Solid: accumulation | Dashed: depletion"
        else:
            legend_text = "Blue: accumulation | Red: depletion"
        ax.text(0.98, 0.02, legend_text,
                transform=ax.transAxes, fontsize=6, va='bottom', ha='right',
                color='gray')

        plt.tight_layout(pad=0.5)
        fig.savefig(pdf_path, format='pdf', bbox_inches='tight')
        plt.close(fig)

        pdf_paths.append(pdf_path)

    return pdf_paths


def combine_pdfs(pdf_paths: list[Path], output_path: Path):
    """Combine multiple PDFs into a single file."""
    from PyPDF2 import PdfMerger

    merger = PdfMerger()
    for pdf in pdf_paths:
        merger.append(str(pdf))
    merger.write(str(output_path))
    merger.close()


def create_info_page(
    molecule_name: str,
    formula: str,
    settings: ContourSettings,
    output_path: Path,
):
    """Create an info page with molecule details and instructions."""
    fig, ax = plt.subplots(figsize=settings.page_size)
    fig.patch.set_facecolor('white')
    ax.axis('off')

    text = f"""
Dorothy - Crystallography Teaching Tool

Molecule: {molecule_name}
Formula: {formula}
Number of slices: {settings.n_slices}

Instructions:
1. Print each slice on transparent film
2. Apply film to 2mm transparent PVC sheets
3. Stack slices in order (1 at bottom, {settings.n_slices} at top)
4. Align using corner marks
5. View the 3D electron density!

Promolecule density shows atomic positions.
Deformation density shows chemical bonding:
  - Solid/Blue: electron accumulation (bonds, lone pairs)
  - Dashed/Red: electron depletion

Inspired by Dorothy Hodgkin's pioneering crystallography work.
    """

    ax.text(0.5, 0.5, text.strip(), transform=ax.transAxes,
            fontsize=10, va='center', ha='center',
            family='monospace')

    fig.savefig(output_path, format='pdf', bbox_inches='tight')
    plt.close(fig)
