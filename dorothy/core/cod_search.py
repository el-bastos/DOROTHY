"""
COD (Crystallography Open Database) search module.

Provides search functionality with fallback to local examples when COD is unavailable.
"""

import requests
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

    def __init__(self, examples_dir: Optional[str | Path] = None):
        self._online = None
        self._local_molecules: list[MoleculeStructure] = []
        self._local_results: list[MoleculeResult] = []

        # Load local examples
        if examples_dir is None:
            # Default: look for examples folder relative to package
            examples_dir = Path(__file__).parent.parent.parent / "examples"

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
                f"{self.BASE_URL}/search.php",
                timeout=5
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

    def _search_cod(self, query: str) -> list[MoleculeResult]:
        """Search COD database online."""
        params = {
            "text": query,
            "format": "json",
        }

        response = requests.get(
            f"{self.BASE_URL}/result.php",
            params=params,
            timeout=self.TIMEOUT
        )
        response.raise_for_status()

        data = response.json()
        results = []

        for entry in data.get("results", [])[:20]:  # Limit to 20 results
            results.append(MoleculeResult(
                cod_id=str(entry.get("file", "")),
                name=entry.get("commonname", entry.get("chemname", "Unknown")),
                formula=entry.get("formula", ""),
                atom_count=entry.get("nel", 0),
                space_group=entry.get("sg", ""),
                is_local=False,
            ))

        return results

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
            # Download from COD
            cif_content = self.download_cif(result.cod_id)
            if cif_content:
                from dorothy.core.cif_parser import parse_cif
                # Would need to save to temp file or parse from string
                # For now, return None for COD results
                pass
        return None

    def get_cif_url(self, cod_id: str) -> str:
        """Get the URL to download a CIF file from COD."""
        return f"{self.BASE_URL}/{cod_id}.cif"

    def download_cif(self, cod_id: str) -> Optional[str]:
        """Download CIF file content from COD."""
        try:
            url = self.get_cif_url(cod_id)
            response = requests.get(url, timeout=self.TIMEOUT)
            response.raise_for_status()
            return response.text
        except (requests.RequestException, requests.Timeout):
            return None
