"""
COD (Crystallography Open Database) search module.

Provides search functionality with fallback to local examples when COD is unavailable.
"""

import warnings
import requests

# Suppress SSL warnings for COD (their certificate has hostname mismatch)
from urllib3.exceptions import InsecureRequestWarning
warnings.filterwarnings('ignore', category=InsecureRequestWarning)
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from dorothy.core.cif_parser import MoleculeStructure, load_example_molecules


@dataclass
class MoleculeResult:
    """Search result - either from COD or local examples."""
    cod_id: str
    name: str
    formula: str
    atom_count: int
    space_group: str = ""
    is_local: bool = False  # True if from local examples
    cif_path: str = ""  # Path to local CIF file (if local)


class CODSearch:
    """Search interface for molecules - COD online + local examples."""

    BASE_URL = "https://www.crystallography.net/cod"
    TIMEOUT = 10  # seconds
    # COD server filters by User-Agent, so we need to use a browser-like one
    HEADERS = {
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) Dorothy/0.3',
        'Accept': '*/*',
    }

    def __init__(self, examples_dir: Optional[str | Path] = None):
        self._online = None
        self._local_molecules: list[MoleculeStructure] = []
        self._local_results: list[MoleculeResult] = []

        # Load local examples
        if examples_dir is None:
            # Default: look for examples folder relative to package
            from dorothy import _base_dir
            examples_dir = _base_dir() / "examples"

        self._load_local_examples(examples_dir)

    def _load_local_examples(self, examples_dir: Path):
        """Load molecules from local examples directory."""
        examples_dir = Path(examples_dir)
        if not examples_dir.exists():
            return

        self._local_molecules = load_example_molecules(examples_dir)
        self._local_results = []

        for mol in self._local_molecules:
            # Extract ID from filename or use name
            cif_path = Path(mol.source_file)
            cod_id = cif_path.stem  # filename without extension

            self._local_results.append(MoleculeResult(
                cod_id=cod_id,
                name=mol.name,
                formula=mol.formula,
                atom_count=mol.atom_count,
                space_group=mol.space_group,
                is_local=True,
                cif_path=mol.source_file,
            ))

    def is_online(self) -> bool:
        """Check if COD is reachable."""
        try:
            response = requests.head(
                f"{self.BASE_URL}/search.html",
                headers=self.HEADERS,
                timeout=5,
                verify=False,  # COD has SSL certificate mismatch
            )
            self._online = response.status_code == 200
        except (requests.RequestException, requests.Timeout):
            self._online = False
        return self._online

    def search(self, query: str) -> tuple[list[MoleculeResult], bool]:
        """
        Search for molecules by name or formula.

        Returns:
            Tuple of (results list, is_online)
            Searches local examples first, then COD if online.
        """
        query_lower = query.lower().strip()
        results = []

        # Always search local examples first
        local_matches = self._search_local(query_lower)
        results.extend(local_matches)

        # Try COD online search
        is_online = False
        try:
            cod_results = self._search_cod(query)
            if cod_results:
                results.extend(cod_results)
                is_online = True
        except (requests.RequestException, requests.Timeout):
            pass

        return results, is_online

    def _search_local(self, query: str) -> list[MoleculeResult]:
        """Search local example molecules."""
        matches = []
        for mol in self._local_results:
            if (query in mol.name.lower() or
                query in mol.formula.lower() or
                query in mol.cod_id.lower()):
                matches.append(mol)
        return matches

    def _estimate_atom_count(self, formula: str) -> int:
        """
        Estimate atom count from chemical formula.

        Parses formulas like "C9 H8 O4" to sum atom counts.
        """
        import re
        total = 0
        # Match element symbol followed by optional number
        # e.g., "C9" -> C with 9, "H" -> H with 1
        for match in re.finditer(r'([A-Z][a-z]?)\s*(\d*)', formula):
            element = match.group(1)
            count_str = match.group(2)
            count = int(count_str) if count_str else 1
            total += count
        return total

    def _search_cod(self, query: str, retries: int = 2) -> list[MoleculeResult]:
        """Search COD database online with retry logic."""
        import time

        last_error = None
        for attempt in range(retries + 1):
            try:
                # COD API requires POST with text1 parameter
                data = {
                    "text1": query,
                    "format": "json",
                }

                response = requests.post(
                    f"{self.BASE_URL}/result.php",
                    data=data,
                    headers=self.HEADERS,
                    timeout=self.TIMEOUT,
                    verify=False,  # COD has SSL certificate mismatch
                )
                response.raise_for_status()

                # COD returns a JSON array directly, not an object with "results" key
                result_data = response.json()
                results = []

                # Handle both array (current API) and object with "results" key (legacy)
                entries = result_data if isinstance(result_data, list) else result_data.get("results", [])

                for entry in entries[:20]:  # Limit to 20 results
                    # Get name: prefer commonname, fall back to chemname
                    name = entry.get("commonname") or entry.get("chemname") or "Unknown"

                    # COD doesn't provide atom count directly
                    # nel = number of element types, Z = molecules per unit cell
                    # We estimate from formula if possible, otherwise use 0
                    formula = entry.get("formula", "").strip("- ")
                    atom_count = self._estimate_atom_count(formula)

                    results.append(MoleculeResult(
                        cod_id=str(entry.get("file", "")),
                        name=name,
                        formula=formula,
                        atom_count=atom_count,
                        space_group=entry.get("sg", ""),
                        is_local=False,
                    ))

                return results

            except (requests.RequestException, requests.Timeout, ValueError) as e:
                last_error = e
                if attempt < retries:
                    time.sleep(0.5 * (attempt + 1))  # Backoff: 0.5s, 1s
                    continue
                raise last_error

    def get_structure(self, result: MoleculeResult) -> Optional[MoleculeStructure]:
        """
        Get full structure for a search result.

        For local results, parses the CIF file.
        For COD results, downloads and parses.
        """
        if result.is_local and result.cif_path:
            # Find in cached local molecules
            for mol in self._local_molecules:
                if mol.source_file == result.cif_path:
                    return mol
        else:
            # Download from COD and parse
            cif_content = self.download_cif(result.cod_id)
            if cif_content:
                from dorothy.core.cif_parser import parse_cif_string
                structure = parse_cif_string(
                    cif_content,
                    name=result.name,
                    source=f"COD:{result.cod_id}"
                )
                return structure
        return None

    def get_cif_url(self, cod_id: str) -> str:
        """Get the URL to download a CIF file from COD."""
        return f"{self.BASE_URL}/{cod_id}.cif"

    def download_cif(self, cod_id: str, retries: int = 2) -> Optional[str]:
        """Download CIF file content from COD with retry logic."""
        import time

        url = self.get_cif_url(cod_id)
        for attempt in range(retries + 1):
            try:
                response = requests.get(
                    url,
                    headers=self.HEADERS,
                    timeout=self.TIMEOUT,
                    verify=False,  # COD has SSL certificate mismatch
                )
                response.raise_for_status()
                return response.text
            except (requests.RequestException, requests.Timeout):
                if attempt < retries:
                    time.sleep(0.5 * (attempt + 1))  # Backoff: 0.5s, 1s
                    continue
                return None
        return None
