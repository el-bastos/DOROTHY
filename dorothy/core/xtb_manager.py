"""
xTB download and execution manager.

Handles automatic download of xTB binaries and running calculations.
"""

import platform
import subprocess
import tempfile
import shutil
import requests
import tarfile
import zipfile
from pathlib import Path
from typing import Optional, Callable

from dorothy.core.cif_parser import MoleculeStructure


# xTB release info
XTB_VERSION = "6.7.1"
XTB_RELEASES = {
    "Darwin": {
        "arm64": f"https://github.com/grimme-lab/xtb/releases/download/v{XTB_VERSION}/xtb-{XTB_VERSION}-macos-arm64.tar.xz",
        "x86_64": f"https://github.com/grimme-lab/xtb/releases/download/v{XTB_VERSION}/xtb-{XTB_VERSION}-macos-x86_64.tar.xz",
    },
    "Linux": {
        "x86_64": f"https://github.com/grimme-lab/xtb/releases/download/v{XTB_VERSION}/xtb-{XTB_VERSION}-linux-x86_64.tar.xz",
    },
    "Windows": {
        "AMD64": f"https://github.com/grimme-lab/xtb/releases/download/v{XTB_VERSION}/xtb-{XTB_VERSION}-windows-x86_64.zip",
    },
}


def get_xtb_dir() -> Path:
    """Get the directory where xTB should be stored."""
    # Use app data directory
    if platform.system() == "Darwin":
        base = Path.home() / "Library" / "Application Support" / "Dorothy"
    elif platform.system() == "Windows":
        base = Path.home() / "AppData" / "Local" / "Dorothy"
    else:
        base = Path.home() / ".local" / "share" / "dorothy"

    base.mkdir(parents=True, exist_ok=True)
    return base / "xtb"


def get_xtb_executable() -> Optional[Path]:
    """Get path to xTB executable if installed."""
    xtb_dir = get_xtb_dir()

    if platform.system() == "Windows":
        exe = xtb_dir / "bin" / "xtb.exe"
    else:
        exe = xtb_dir / "bin" / "xtb"

    if exe.exists():
        return exe
    return None


def is_xtb_installed() -> bool:
    """Check if xTB is installed."""
    return get_xtb_executable() is not None


def get_download_url() -> Optional[str]:
    """Get the xTB download URL for current platform."""
    system = platform.system()
    machine = platform.machine()

    if system not in XTB_RELEASES:
        return None

    arch_map = XTB_RELEASES[system]
    return arch_map.get(machine)


def download_xtb(progress_callback: Optional[Callable[[int, int], None]] = None) -> bool:
    """
    Download and install xTB.

    Args:
        progress_callback: Optional callback(downloaded, total) for progress updates

    Returns:
        True if successful, False otherwise
    """
    url = get_download_url()
    if not url:
        return False

    xtb_dir = get_xtb_dir()

    try:
        # Download to temp file
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))

        with tempfile.NamedTemporaryFile(delete=False, suffix='.tar.xz') as tmp:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                tmp.write(chunk)
                downloaded += len(chunk)
                if progress_callback:
                    progress_callback(downloaded, total_size)
            tmp_path = Path(tmp.name)

        # Extract
        if xtb_dir.exists():
            shutil.rmtree(xtb_dir)
        xtb_dir.mkdir(parents=True)

        if url.endswith('.zip'):
            with zipfile.ZipFile(tmp_path, 'r') as zf:
                zf.extractall(xtb_dir.parent)
        else:
            # tar.xz
            with tarfile.open(tmp_path, 'r:xz') as tf:
                tf.extractall(xtb_dir.parent)

        # Rename extracted folder to 'xtb'
        extracted = list(xtb_dir.parent.glob(f'xtb-{XTB_VERSION}*'))
        if extracted:
            if xtb_dir.exists():
                shutil.rmtree(xtb_dir)
            extracted[0].rename(xtb_dir)

        # Make executable on Unix
        if platform.system() != "Windows":
            exe = xtb_dir / "bin" / "xtb"
            if exe.exists():
                exe.chmod(0o755)

        # Cleanup
        tmp_path.unlink()

        return is_xtb_installed()

    except Exception as e:
        print(f"Failed to download xTB: {e}")
        return False


def write_xyz_file(structure: MoleculeStructure, filepath: Path):
    """Write structure to XYZ format for xTB input."""
    coords = structure.get_cartesian_coords()
    symbols = structure.get_symbols()

    with open(filepath, 'w') as f:
        f.write(f"{len(symbols)}\n")
        f.write(f"{structure.name}\n")
        for sym, (x, y, z) in zip(symbols, coords):
            f.write(f"{sym:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")


def run_xtb_density(
    structure: MoleculeStructure,
    output_dir: Path,
    progress_callback: Optional[Callable[[str], None]] = None
) -> Optional[Path]:
    """
    Run xTB calculation to generate electron density cube file.

    Args:
        structure: Molecule structure
        output_dir: Directory for output files
        progress_callback: Optional callback for status updates

    Returns:
        Path to density cube file, or None if failed
    """
    exe = get_xtb_executable()
    if not exe:
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write input XYZ file
    xyz_file = output_dir / "molecule.xyz"
    write_xyz_file(structure, xyz_file)

    if progress_callback:
        progress_callback("Running xTB calculation...")

    try:
        # Run xTB with density output
        # --cube generates electron density cube file
        result = subprocess.run(
            [str(exe), str(xyz_file), "--gfn", "2", "--cube"],
            cwd=output_dir,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
        )

        if result.returncode != 0:
            print(f"xTB error: {result.stderr}")
            return None

        # xTB outputs density.cub
        cube_file = output_dir / "density.cub"
        if cube_file.exists():
            return cube_file

        # Alternative name
        cube_file = output_dir / "xtb_density.cub"
        if cube_file.exists():
            return cube_file

        return None

    except subprocess.TimeoutExpired:
        print("xTB calculation timed out")
        return None
    except Exception as e:
        print(f"xTB error: {e}")
        return None
