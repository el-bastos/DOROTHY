"""
Main generation pipeline for Dorothy.

Orchestrates the full workflow:
1. Check/download xTB
2. Run density calculation
3. Generate contour slices
4. Export PDFs
"""

import shutil
import tempfile
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Callable

from dorothy.core.cif_parser import MoleculeStructure
from dorothy.core.xtb_manager import is_xtb_installed, run_xtb_density
from dorothy.core.density import (
    parse_cube_file,
    calculate_deformation_density,
    create_density_cube_from_structure,
    DensityCube,
)
from dorothy.core.contours import (
    ContourSettings,
    generate_promolecule_contours,
    generate_deformation_contours,
    create_info_page,
)


@dataclass
class GenerationSettings:
    """Settings for the generation pipeline."""
    resolution: str = "medium"  # coarse, medium, fine
    n_slices: int = 15
    generate_promolecule: bool = True
    generate_deformation: bool = True
    color_mode: str = "bw"  # bw, color


@dataclass
class GenerationResult:
    """Result of the generation pipeline."""
    success: bool
    output_dir: Optional[Path] = None
    promolecule_pdfs: list[Path] = None
    deformation_pdfs: list[Path] = None
    error_message: str = ""
    used_xtb: bool = False

    def __post_init__(self):
        if self.promolecule_pdfs is None:
            self.promolecule_pdfs = []
        if self.deformation_pdfs is None:
            self.deformation_pdfs = []


class GenerationPipeline:
    """
    Main pipeline for generating electron density contour slices.
    """

    def __init__(
        self,
        progress_callback: Optional[Callable[[str, int], None]] = None
    ):
        """
        Args:
            progress_callback: Callback(message, percent) for progress updates
        """
        self.progress_callback = progress_callback

    def _report_progress(self, message: str, percent: int):
        """Report progress to callback."""
        if self.progress_callback:
            self.progress_callback(message, percent)

    def generate(
        self,
        structure: MoleculeStructure,
        output_dir: Path,
        settings: GenerationSettings,
    ) -> GenerationResult:
        """
        Run the full generation pipeline.

        Args:
            structure: Molecule structure
            output_dir: Directory for output files
            settings: Generation settings

        Returns:
            GenerationResult with output file paths
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        result = GenerationResult(success=False, output_dir=output_dir)

        try:
            # Step 1: Get molecular density
            self._report_progress("Calculating electron density...", 10)

            molecular_density: Optional[DensityCube] = None
            promolecule_density: Optional[DensityCube] = None

            if is_xtb_installed():
                # Use xTB for accurate molecular density
                self._report_progress("Running xTB calculation...", 20)

                with tempfile.TemporaryDirectory() as tmp_dir:
                    cube_path = run_xtb_density(
                        structure,
                        Path(tmp_dir),
                        lambda msg: self._report_progress(msg, 30)
                    )

                    if cube_path:
                        molecular_density = parse_cube_file(cube_path)
                        result.used_xtb = True

            if molecular_density is None:
                # Fall back to promolecule-only mode
                self._report_progress("Generating promolecule density...", 30)
                promolecule_density = create_density_cube_from_structure(
                    structure, settings.resolution
                )
            else:
                # Create promolecule on same grid as molecular
                self._report_progress("Calculating promolecule density...", 40)
                promolecule_density = create_density_cube_from_structure(
                    structure, settings.resolution
                )

            # Step 2: Calculate deformation density if we have molecular density
            deformation_density: Optional[DensityCube] = None
            if molecular_density and settings.generate_deformation:
                self._report_progress("Computing deformation density...", 50)
                deformation_density = calculate_deformation_density(
                    molecular_density, structure
                )

            # Step 3: Generate contour PDFs
            contour_settings = ContourSettings(
                n_slices=settings.n_slices,
                color_mode=settings.color_mode,
            )

            # Promolecule contours
            if settings.generate_promolecule and promolecule_density:
                self._report_progress("Generating promolecule slices...", 60)
                promol_dir = output_dir / "promolecule"

                def promol_progress(curr, total, msg):
                    pct = 60 + int(20 * curr / total)
                    self._report_progress(msg, pct)

                result.promolecule_pdfs = generate_promolecule_contours(
                    promolecule_density,
                    promol_dir,
                    contour_settings,
                    promol_progress,
                )

            # Deformation contours
            if settings.generate_deformation and deformation_density:
                self._report_progress("Generating deformation slices...", 80)
                deform_dir = output_dir / "deformation"

                def deform_progress(curr, total, msg):
                    pct = 80 + int(15 * curr / total)
                    self._report_progress(msg, pct)

                result.deformation_pdfs = generate_deformation_contours(
                    deformation_density,
                    deform_dir,
                    contour_settings,
                    deform_progress,
                )

            # Step 4: Create info page
            self._report_progress("Creating info page...", 95)
            info_path = output_dir / "info.pdf"
            create_info_page(
                structure.name,
                structure.formula,
                contour_settings,
                info_path,
            )

            self._report_progress("Complete!", 100)
            result.success = True

        except Exception as e:
            result.error_message = str(e)
            self._report_progress(f"Error: {e}", 0)

        return result
