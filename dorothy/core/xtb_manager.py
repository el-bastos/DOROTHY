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
from typing import Optional, Callable, TYPE_CHECKING

import numpy as np

from dorothy.core.cif_parser import MoleculeStructure

if TYPE_CHECKING:
    from dorothy.core.selection import PlaneDefinition


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
    # First check our managed installation
    xtb_dir = get_xtb_dir()

    if platform.system() == "Windows":
        exe = xtb_dir / "bin" / "xtb.exe"
    else:
        exe = xtb_dir / "bin" / "xtb"

    if exe.exists():
        return exe

    # Check if xTB is in PATH (homebrew, conda, etc.)
    xtb_in_path = shutil.which("xtb")
    if xtb_in_path:
        return Path(xtb_in_path)

    # Check common installation locations (PATH may not be set in GUI apps)
    common_locations = []
    if platform.system() == "Darwin":
        # Homebrew locations
        common_locations = [
            Path("/opt/homebrew/bin/xtb"),  # Apple Silicon
            Path("/usr/local/bin/xtb"),      # Intel Mac
            Path.home() / "miniconda3" / "bin" / "xtb",
            Path.home() / "anaconda3" / "bin" / "xtb",
            Path.home() / "miniforge3" / "bin" / "xtb",
        ]
    elif platform.system() == "Linux":
        common_locations = [
            Path("/usr/local/bin/xtb"),
            Path("/usr/bin/xtb"),
            Path.home() / "miniconda3" / "bin" / "xtb",
            Path.home() / "anaconda3" / "bin" / "xtb",
            Path.home() / "miniforge3" / "bin" / "xtb",
        ]

    for loc in common_locations:
        if loc.exists():
            return loc

    return None


def is_xtb_installed() -> bool:
    """Check if xTB is installed."""
    return get_xtb_executable() is not None


def get_xtb_install_instructions() -> str:
    """Get platform-specific installation instructions."""
    system = platform.system()

    if system == "Darwin":
        return (
            "xTB can be installed via Homebrew:\n\n"
            "  brew tap grimme-lab/qc\n"
            "  brew install xtb\n\n"
            "Or via conda:\n\n"
            "  conda install -c conda-forge xtb"
        )
    elif system == "Linux":
        return (
            "xTB can be installed via conda:\n\n"
            "  conda install -c conda-forge xtb\n\n"
            "Or download from GitHub releases."
        )
    elif system == "Windows":
        return (
            "xTB can be installed via conda:\n\n"
            "  conda install -c conda-forge xtb\n\n"
            "Or download from GitHub releases."
        )
    else:
        return "Please install xTB from https://github.com/grimme-lab/xtb"


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


def write_xyz_file(
    structure: MoleculeStructure,
    filepath: Path,
    align_to_principal_axes: bool = True,
    plane_definition: Optional["PlaneDefinition"] = None
):
    """Write structure to XYZ format for xTB input.

    Args:
        structure: Molecule structure
        filepath: Path to write XYZ file
        align_to_principal_axes: If True and no plane_definition, align to principal axes
        plane_definition: If provided, rotate coordinates so the user-selected plane
                         becomes horizontal (XY plane). This takes precedence over
                         align_to_principal_axes.
    """
    if plane_definition is not None:
        # Use custom plane rotation: rotate so plane normal -> Z axis
        # IMPORTANT: The plane_definition was calculated from principal-axes-aligned
        # coordinates, so we must start from those same aligned coordinates
        coords = structure.get_cartesian_coords(align_to_principal_axes=True)
        center = plane_definition.center
        rotation = plane_definition.rotation_matrix
        # Center on plane center, then rotate
        coords_centered = coords - center
        coords = coords_centered @ rotation.T
    else:
        coords = structure.get_cartesian_coords(align_to_principal_axes=align_to_principal_axes)

    symbols = structure.get_symbols()

    with open(filepath, 'w') as f:
        f.write(f"{len(symbols)}\n")
        f.write(f"{structure.name}\n")
        for sym, (x, y, z) in zip(symbols, coords):
            f.write(f"{sym:2s} {x:15.8f} {y:15.8f} {z:15.8f}\n")


def run_xtb_density(
    structure: MoleculeStructure,
    output_dir: Path,
    progress_callback: Optional[Callable[[str], None]] = None,
    plane_definition: Optional["PlaneDefinition"] = None,
) -> Optional[Path]:
    """
    Run xTB calculation to generate electron density cube file.

    Args:
        structure: Molecule structure
        output_dir: Directory for output files
        progress_callback: Optional callback for status updates
        plane_definition: If provided, rotate molecule so the user-selected plane
                         becomes horizontal before calculation. This ensures
                         Z-slices cut parallel to the user's chosen plane.

    Returns:
        Path to density cube file, or None if failed
    """
    exe = get_xtb_executable()
    if not exe:
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write input XYZ file (with optional plane rotation)
    xyz_file = output_dir / "molecule.xyz"
    write_xyz_file(structure, xyz_file, plane_definition=plane_definition)

    # Write xTB input file to request density cube output
    input_file = output_dir / "xtb.inp"
    input_file.write_text("$write\n   density=true\n$end\n")

    if progress_callback:
        progress_callback("Running xTB calculation...")

    try:
        # Run xTB with density output via input file
        result = subprocess.run(
            [str(exe), str(xyz_file), "--gfn", "2", "--input", str(input_file)],
            cwd=output_dir,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
        )

        if result.returncode != 0:
            print(f"xTB error: {result.stderr}")
            return None

        # xTB outputs density.cub when density=true is set
        density_file = output_dir / "density.cub"
        if not density_file.exists():
            # Check for other possible cube file names
            for pattern in ["*density*.cub", "*density*.cube"]:
                cubes = list(output_dir.glob(pattern))
                if cubes:
                    density_file = cubes[0]
                    break
            else:
                density_file = None

        return density_file

    except subprocess.TimeoutExpired:
        print("xTB calculation timed out")
        return None
    except Exception as e:
        print(f"xTB error: {e}")
        return None
