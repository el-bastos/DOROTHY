"""
CIF (Crystallographic Information File) parser.

Extracts molecular structure data from CIF files.
Handles symmetry expansion to reconstruct complete molecules from asymmetric units.
"""

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

logger = logging.getLogger(__name__)

# Covalent radii (Å) for bond detection during molecule extraction.
# Source: Cordero et al., Dalton Trans., 2008, 2832-2838
_COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07,
    'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76,
    'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.39,
    'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
    'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20,
    'Rb': 2.20, 'Sr': 1.95, 'Zr': 1.75, 'Mo': 1.54, 'Ru': 1.46,
    'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42,
    'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Cs': 2.44,
    'Ba': 2.15, 'La': 2.07, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
    'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48,
}
_DEFAULT_COVALENT_RADIUS = 1.5


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


def _parse_symop_component(expr: str) -> dict:
    """Parse one component of a symmetry operation (e.g., '-x+1/2').

    Returns dict with 'x', 'y', 'z' coefficients and 'c' constant.
    Example: '-x+1/2' → {'x': -1.0, 'y': 0.0, 'z': 0.0, 'c': 0.5}
    """
    result = {'x': 0.0, 'y': 0.0, 'z': 0.0, 'c': 0.0}
    expr = expr.strip().replace(' ', '')

    tokens = re.findall(r'[+-]?[^+-]+', expr)

    for token in tokens:
        token = token.strip()
        if not token:
            continue

        found_var = False
        for var in ('x', 'y', 'z'):
            if var in token:
                coeff_str = token.replace(var, '')
                if coeff_str in ('', '+'):
                    result[var] = 1.0
                elif coeff_str == '-':
                    result[var] = -1.0
                else:
                    result[var] = float(coeff_str)
                found_var = True
                break

        if not found_var:
            # Constant term (number or fraction)
            if '/' in token:
                num, den = token.split('/')
                result['c'] += float(num) / float(den)
            else:
                result['c'] += float(token)

    return result


def _parse_symmetry_ops(content: str) -> list[np.ndarray]:
    """Parse symmetry operations from CIF content.

    Looks for _space_group_symop_operation_xyz or _symmetry_equiv_pos_as_xyz.

    Returns list of 3×4 numpy arrays: [[Rx, Ry, Rz, T], ...] per operation.
    """
    for tag in ('_space_group_symop_operation_xyz', '_symmetry_equiv_pos_as_xyz'):
        if tag not in content:
            continue

        tag_pos = content.index(tag)
        remaining = content[tag_pos:]

        ops = []
        for line in remaining.split('\n')[1:]:  # skip the header line
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('_') or stripped.startswith('loop_'):
                break

            # Extract quoted operation string, or use whole line
            op_match = re.search(r"'([^']+)'|\"([^\"]+)\"", stripped)
            if op_match:
                op_str = op_match.group(1) or op_match.group(2)
            else:
                op_str = stripped
                # Skip leading numeric ID if present (e.g., "1 x, y, z")
                op_str = re.sub(r'^\d+\s+', '', op_str)

            components = [c.strip() for c in op_str.split(',')]
            if len(components) != 3:
                continue

            # Verify it looks like a symmetry operation
            if not any(v in ''.join(components) for v in ('x', 'y', 'z')):
                continue

            try:
                parsed = [_parse_symop_component(c) for c in components]
                matrix = np.array([
                    [p['x'], p['y'], p['z'], p['c']] for p in parsed
                ])
                ops.append(matrix)
            except (ValueError, KeyError):
                continue

        if ops:
            logger.info("Parsed %d symmetry operations", len(ops))
            return ops

    return []


def _apply_symmetry_ops(
    atoms: list['Atom'], sym_ops: list[np.ndarray], tolerance: float = 0.02
) -> list['Atom']:
    """Apply symmetry operations to expand asymmetric unit to full unit cell.

    Args:
        atoms: Asymmetric unit atoms
        sym_ops: List of 3×4 affine transformation matrices
        tolerance: Fractional-coordinate distance for duplicate removal

    Returns:
        All unique atoms in the unit cell
    """
    if not sym_ops:
        return list(atoms)

    all_atoms: list['Atom'] = []
    all_frac: list[np.ndarray] = []

    for atom in atoms:
        frac = np.array([atom.x, atom.y, atom.z])

        for op in sym_ops:
            new_frac = op[:, :3] @ frac + op[:, 3]
            new_frac = new_frac % 1.0

            # Deduplicate (same position within tolerance)
            is_dup = False
            for existing_frac in all_frac:
                diff = existing_frac - new_frac
                diff = diff - np.round(diff)  # periodic distance
                if np.linalg.norm(diff) < tolerance:
                    is_dup = True
                    break

            if not is_dup:
                all_atoms.append(Atom(
                    label=atom.label,
                    symbol=atom.symbol,
                    x=new_frac[0],
                    y=new_frac[1],
                    z=new_frac[2],
                ))
                all_frac.append(new_frac)

    logger.info(
        "Symmetry expansion: %d asymmetric → %d unit cell atoms",
        len(atoms), len(all_atoms),
    )
    return all_atoms


def _extract_molecule(
    atoms: list['Atom'], cell: 'CellParameters',
    bond_tolerance: float = 1.3, z_units: int = 0,
) -> list['Atom']:
    """Extract one formula unit from unit cell atoms using bond connectivity.

    Detects bonds via covalent radii, finds connected components (considering
    periodic boundary conditions), picks the largest fragment, then merges
    nearby counterions/solvate until the expected formula-unit size is reached.

    Args:
        atoms: All atoms in the unit cell
        cell: Unit cell parameters
        bond_tolerance: Multiplier for sum of covalent radii
        z_units: Number of formula units per unit cell (from _cell_formula_units_Z).
                 If > 0, used to determine expected atom count per formula unit.

    Returns:
        Atoms of one complete formula unit
    """
    n = len(atoms)
    if n <= 1:
        return atoms

    matrix = cell.to_cartesian_matrix()
    frac_coords = np.array([[a.x, a.y, a.z] for a in atoms])
    cart_coords = frac_coords @ matrix.T

    radii = [_COVALENT_RADII.get(a.symbol, _DEFAULT_COVALENT_RADIUS) for a in atoms]

    # Build adjacency considering periodic images (27 cells)
    adjacency: list[set[int]] = [set() for _ in range(n)]
    bond_offsets: dict[tuple[int, int], np.ndarray] = {}

    for i in range(n):
        for j in range(i + 1, n):
            max_dist = (radii[i] + radii[j]) * bond_tolerance
            best_dist = float('inf')
            best_offset = None

            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        offset = np.array([dx, dy, dz], dtype=float)
                        shifted_cart = (frac_coords[j] + offset) @ matrix.T
                        dist = np.linalg.norm(cart_coords[i] - shifted_cart)
                        if 0.4 < dist < max_dist and dist < best_dist:
                            best_dist = dist
                            best_offset = offset

            if best_offset is not None:
                adjacency[i].add(j)
                adjacency[j].add(i)
                bond_offsets[(i, j)] = best_offset
                bond_offsets[(j, i)] = -best_offset

    # BFS to find connected components, tracking offsets for unwrapping
    visited = [False] * n
    components: list[tuple[list[int], dict[int, np.ndarray]]] = []

    for start in range(n):
        if visited[start]:
            continue
        component: list[int] = []
        offsets: dict[int, np.ndarray] = {start: np.zeros(3)}
        queue = [start]
        visited[start] = True

        while queue:
            current = queue.pop(0)
            component.append(current)

            for neighbor in adjacency[current]:
                if not visited[neighbor]:
                    visited[neighbor] = True
                    queue.append(neighbor)
                    if (current, neighbor) in bond_offsets:
                        offsets[neighbor] = (
                            offsets[current] + bond_offsets[(current, neighbor)]
                        )
                    else:
                        offsets[neighbor] = offsets[current].copy()

        components.append((component, offsets))

    # Pick the largest connected component as the starting fragment.
    best_idx = max(range(len(components)), key=lambda i: len(components[i][0]))
    best_comp, best_offsets = components[best_idx]

    molecule_indices = list(best_comp)
    molecule_offsets = dict(best_offsets)

    # Determine expected formula-unit size from Z.
    # If Z is known, merge nearby fragments (counterions, solvate) until
    # we reach the expected atom count.
    target_size = (n // z_units) if z_units > 0 else 0

    if target_size > len(molecule_indices):
        remaining = [
            (i, c) for i, c in enumerate(components) if i != best_idx
        ]

        # Compute Cartesian coords of the growing molecule
        mol_frac = np.array([
            frac_coords[idx] + molecule_offsets[idx]
            for idx in molecule_indices
        ])
        mol_cart = mol_frac @ matrix.T

        while remaining and len(molecule_indices) < target_size:
            # Find the closest remaining fragment (checking 27 periodic images)
            closest_i = None
            closest_dist = float('inf')
            closest_shift = np.zeros(3)

            for ri, (_, (comp, comp_offs)) in enumerate(remaining):
                comp_frac = np.array([
                    frac_coords[idx] + comp_offs[idx] for idx in comp
                ])
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for dz in (-1, 0, 1):
                            shift = np.array([dx, dy, dz], dtype=float)
                            shifted_cart = (comp_frac + shift) @ matrix.T
                            diffs = mol_cart[:, None, :] - shifted_cart[None, :, :]
                            dists = np.linalg.norm(diffs, axis=2)
                            min_d = dists.min()
                            if min_d < closest_dist:
                                closest_dist = min_d
                                closest_i = ri
                                closest_shift = shift

            if closest_i is None:
                break

            # Merge the closest fragment
            _, (comp, comp_offs) = remaining.pop(closest_i)
            for idx in comp:
                molecule_offsets[idx] = comp_offs[idx] + closest_shift
            molecule_indices.extend(comp)

            # Update mol_cart for next iteration
            mol_frac = np.array([
                frac_coords[idx] + molecule_offsets[idx]
                for idx in molecule_indices
            ])
            mol_cart = mol_frac @ matrix.T

    logger.info(
        "Molecule extraction: %d unit-cell atoms → %d fragments, "
        "formula unit has %d atoms (Z=%s, target=%s)",
        n, len(components), len(molecule_indices),
        z_units or '?', target_size or 'largest',
    )

    # Build molecule with unwrapped fractional coordinates
    result = []
    for idx in molecule_indices:
        a = atoms[idx]
        offset = molecule_offsets.get(idx, np.zeros(3))
        result.append(Atom(
            label=a.label,
            symbol=a.symbol,
            x=a.x + offset[0],
            y=a.y + offset[1],
            z=a.z + offset[2],
        ))

    return result


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

    # Space group — try multiple CIF tags (older and newer conventions)
    for sg_tag in (
        r'_symmetry_space_group_name_H-M',
        r'_space_group_name_H-M_alt',
    ):
        sg_match = re.search(rf"{sg_tag}\s+'([^']+)'", content)
        if not sg_match:
            sg_match = re.search(rf'{sg_tag}\s+(\S+)', content)
        if sg_match:
            structure.space_group = sg_match.group(1)
            break

    # Formula units per unit cell (Z) — used for ionic compound extraction
    z_units = 0
    z_match = re.search(r'_cell_formula_units_Z\s+(\d+)', content)
    if z_match:
        z_units = int(z_match.group(1))

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
            occ_idx = headers.index('occupancy') if 'occupancy' in headers else None
            min_col = max(label_idx, symbol_idx, x_idx, y_idx, z_idx)

            atoms_with_occ: list[tuple[Atom, float]] = []
            for row in rows:
                if len(row) > min_col:
                    occ = 1.0
                    if occ_idx is not None and occ_idx < len(row):
                        occ = parse_cif_value(row[occ_idx])
                    atom = Atom(
                        label=row[label_idx],
                        symbol=row[symbol_idx],
                        x=parse_cif_value(row[x_idx]),
                        y=parse_cif_value(row[y_idx]),
                        z=parse_cif_value(row[z_idx]),
                    )
                    atoms_with_occ.append((atom, occ))

            # Resolve disorder: for partial-occupancy atoms, group by element
            # and keep round(sum_of_occupancies) representatives per group.
            # E.g. 8 F atoms at 0.25 occ → sum=2.0 → keep 2 representatives.
            full_occ = [(a, o) for a, o in atoms_with_occ if o >= 0.9]
            partial_by_element: dict[str, list[tuple[Atom, float]]] = {}
            for a, o in atoms_with_occ:
                if o < 0.9:
                    partial_by_element.setdefault(a.symbol, []).append((a, o))

            structure.atoms = [a for a, _ in full_occ]
            for element, group in partial_by_element.items():
                total_occ = sum(o for _, o in group)
                n_keep = max(1, round(total_occ))
                structure.atoms.extend(a for a, _ in group[:n_keep])
                if n_keep < len(group):
                    logger.info(
                        "Disorder resolution: kept %d of %d %s atoms "
                        "(total occupancy %.2f)",
                        n_keep, len(group), element, total_occ,
                    )
        except ValueError:
            # Required columns not found
            pass

    # --- Symmetry expansion ---
    # Apply space group symmetry to reconstruct the full unit cell,
    # then extract one complete molecule via bond connectivity.
    sym_ops = _parse_symmetry_ops(content)
    if sym_ops and structure.atoms:
        asymmetric_count = len(structure.atoms)
        unit_cell_atoms = _apply_symmetry_ops(structure.atoms, sym_ops)
        molecule_atoms = _extract_molecule(
            unit_cell_atoms, structure.cell, z_units=z_units,
        )
        structure.atoms = molecule_atoms
        logger.info(
            "Symmetry: %d asymmetric → %d molecule atoms",
            asymmetric_count, len(molecule_atoms),
        )

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
