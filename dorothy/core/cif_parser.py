"""
CIF (Crystallographic Information File) parser.

Extracts molecular structure data from CIF files.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np


@dataclass
class Atom:
    """Single atom in a structure."""
    label: str
    symbol: str
    x: float  # fractional coordinate
    y: float
    z: float


@dataclass
class CellParameters:
    """Unit cell parameters."""
    a: float
    b: float
    c: float
    alpha: float  # degrees
    beta: float
    gamma: float

    def to_cartesian_matrix(self) -> np.ndarray:
        """Convert cell parameters to Cartesian transformation matrix."""
        # Convert angles to radians
        alpha_rad = np.radians(self.alpha)
        beta_rad = np.radians(self.beta)
        gamma_rad = np.radians(self.gamma)

        # Calculate volume factor
        cos_alpha = np.cos(alpha_rad)
        cos_beta = np.cos(beta_rad)
        cos_gamma = np.cos(gamma_rad)
        sin_gamma = np.sin(gamma_rad)

        # Transformation matrix (fractional -> Cartesian)
        # Standard crystallographic convention
        omega = np.sqrt(1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2
                       + 2 * cos_alpha * cos_beta * cos_gamma)

        matrix = np.array([
            [self.a, self.b * cos_gamma, self.c * cos_beta],
            [0, self.b * sin_gamma, self.c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
            [0, 0, self.c * omega / sin_gamma]
        ])
        return matrix


@dataclass
class MoleculeStructure:
    """Complete molecular structure from CIF."""
    name: str
    formula: str
    cell: CellParameters
    space_group: str
    atoms: list[Atom] = field(default_factory=list)
    source_file: str = ""

    @property
    def atom_count(self) -> int:
        return len(self.atoms)

    def get_cartesian_coords(self, align_to_principal_axes: bool = False) -> np.ndarray:
        """
        Get Cartesian coordinates for all atoms.

        Args:
            align_to_principal_axes: If True, rotate molecule so that:
                - The molecular plane lies in the XY plane
                - The longest axis is along X
                - Z-slicing will cut parallel to the molecular plane
        """
        if not self.atoms:
            return np.array([])

        frac_coords = np.array([[a.x, a.y, a.z] for a in self.atoms])
        matrix = self.cell.to_cartesian_matrix()
        coords = frac_coords @ matrix.T

        if align_to_principal_axes:
            coords = self._align_to_principal_axes(coords)

        return coords

    def _align_to_principal_axes(self, coords: np.ndarray) -> np.ndarray:
        """
        Rotate coordinates to align with principal axes.

        Uses SVD to find the best-fit plane and orients the molecule so:
        - XY plane contains the molecular plane (maximum spread)
        - Z axis is perpendicular to the molecular plane (minimum spread)
        - X axis is along the longest dimension
        """
        # Center the molecule
        center = coords.mean(axis=0)
        centered = coords - center

        # SVD gives principal axes
        # U: left singular vectors (not used)
        # s: singular values (spread along each axis)
        # Vh: right singular vectors (principal directions)
        _, s, vh = np.linalg.svd(centered)

        # vh rows are principal axes in order of decreasing spread
        # vh[0] = direction of maximum spread (longest)
        # vh[1] = direction of medium spread
        # vh[2] = direction of minimum spread (normal to molecular plane)

        # Build rotation matrix: new axes are the principal axes
        # We want: X = longest, Y = medium, Z = shortest (plane normal)
        rotation_matrix = vh  # Each row becomes new axis direction

        # Ensure right-handed coordinate system
        if np.linalg.det(rotation_matrix) < 0:
            rotation_matrix[2] = -rotation_matrix[2]

        # Rotate coordinates
        rotated = centered @ rotation_matrix.T

        # Re-center to positive coordinates (for grid generation)
        rotated = rotated - rotated.min(axis=0)

        return rotated

    def get_symbols(self) -> list[str]:
        """Get list of element symbols."""
        return [a.symbol for a in self.atoms]


def parse_cif_value(value: str) -> float:
    """Parse a CIF numeric value, handling uncertainties like '11.233(3)'."""
    if not value or value == '?':
        return 0.0
    # Remove uncertainty in parentheses
    clean = re.sub(r'\([^)]*\)', '', value)
    try:
        return float(clean)
    except ValueError:
        return 0.0


def parse_cif_string(content: str, name: str = "", source: str = "") -> MoleculeStructure:
    """
    Parse CIF content from a string.

    Args:
        content: CIF file content as string
        name: Optional name for the structure
        source: Optional source identifier (e.g., COD ID)

    Returns:
        MoleculeStructure with all extracted data
    """
    return _parse_cif_content(content, name, source)


def parse_cif(filepath: str | Path) -> MoleculeStructure:
    """
    Parse a CIF file and extract molecular structure.

    Args:
        filepath: Path to CIF file

    Returns:
        MoleculeStructure with all extracted data
    """
    filepath = Path(filepath)
    content = filepath.read_text()
    return _parse_cif_content(content, "", str(filepath))


def _parse_cif_content(content: str, name: str = "", source: str = "") -> MoleculeStructure:
    """
    Internal function to parse CIF content.

    Args:
        content: CIF file content as string
        name: Optional name override
        source: Source file path or identifier

    Returns:
        MoleculeStructure with all extracted data
    """
    # Initialize structure
    structure = MoleculeStructure(
        name=name,
        formula="",
        cell=CellParameters(1, 1, 1, 90, 90, 90),
        space_group="",
        source_file=source
    )

    # Extract name from CIF only if not provided
    if not structure.name:
        name_match = re.search(r'_chemical_name_common\s+(.+)', content)
        if name_match:
            structure.name = name_match.group(1).strip().strip("'\"")

    if not structure.name:
        name_match = re.search(r'_chemical_name_systematic\s*\n;\s*\n(.+?)\n', content)
        if name_match:
            structure.name = name_match.group(1).strip()

    formula_match = re.search(r"_chemical_formula_sum\s+'([^']+)'", content)
    if formula_match:
        structure.formula = formula_match.group(1)

    # Space group can be quoted (e.g., 'P 1 21/c 1') or unquoted
    sg_match = re.search(r"_symmetry_space_group_name_H-M\s+'([^']+)'", content)
    if not sg_match:
        sg_match = re.search(r'_symmetry_space_group_name_H-M\s+(\S+)', content)
    if sg_match:
        structure.space_group = sg_match.group(1)

    # Cell parameters
    for param, attr in [
        ('_cell_length_a', 'a'),
        ('_cell_length_b', 'b'),
        ('_cell_length_c', 'c'),
        ('_cell_angle_alpha', 'alpha'),
        ('_cell_angle_beta', 'beta'),
        ('_cell_angle_gamma', 'gamma'),
    ]:
        match = re.search(rf'{param}\s+(\S+)', content)
        if match:
            setattr(structure.cell, attr, parse_cif_value(match.group(1)))

    # Parse atom sites (this is the tricky part - CIF loop format)
    # We need to find the loop_ block with fractional coordinates specifically
    # and stop before the next loop_ or before lines starting with _
    atom_loop_match = re.search(
        r'loop_\s*\n((?:_atom_site_(?!aniso)\w+\s*\n)+)',
        content,
        re.MULTILINE
    )

    if atom_loop_match:
        # Get column headers
        headers_text = atom_loop_match.group(1)
        headers = re.findall(r'_atom_site_(\w+)', headers_text)

        # Check if this loop has fractional coordinates
        if 'fract_x' not in headers:
            atom_loop_match = None

    if atom_loop_match:
        # Find where the data starts (right after headers)
        data_start = atom_loop_match.end()

        # Extract data rows - stop at next loop_, line starting with _, or empty line followed by _
        data_lines = []
        remaining_content = content[data_start:]
        for line in remaining_content.split('\n'):
            stripped = line.strip()
            # Stop conditions: empty line, line starting with _, or loop_
            if not stripped:
                continue
            if stripped.startswith('_') or stripped.startswith('loop_'):
                break
            data_lines.append(stripped)

        rows = [line.split() for line in data_lines if line]

        # Find column indices
        try:
            label_idx = headers.index('label')
            symbol_idx = headers.index('type_symbol')
            x_idx = headers.index('fract_x')
            y_idx = headers.index('fract_y')
            z_idx = headers.index('fract_z')

            for row in rows:
                if len(row) > max(label_idx, symbol_idx, x_idx, y_idx, z_idx):
                    atom = Atom(
                        label=row[label_idx],
                        symbol=row[symbol_idx],
                        x=parse_cif_value(row[x_idx]),
                        y=parse_cif_value(row[y_idx]),
                        z=parse_cif_value(row[z_idx]),
                    )
                    structure.atoms.append(atom)
        except ValueError:
            # Required columns not found
            pass

    return structure


def load_example_molecules(examples_dir: str | Path) -> list[MoleculeStructure]:
    """
    Load all example molecules from a directory.

    Expects structure: examples_dir/molecule_name/*.cif
    """
    examples_dir = Path(examples_dir)
    molecules = []

    if not examples_dir.exists():
        return molecules

    for mol_dir in examples_dir.iterdir():
        if mol_dir.is_dir() and not mol_dir.name.startswith('.'):
            for cif_file in mol_dir.glob('*.cif'):
                try:
                    structure = parse_cif(cif_file)
                    if not structure.name:
                        structure.name = mol_dir.name.title()
                    molecules.append(structure)
                except Exception as e:
                    print(f"Warning: Failed to parse {cif_file}: {e}")

    return molecules
