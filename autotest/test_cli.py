import platform
import subprocess

from conftest import PROJ_ROOT

bin_path = PROJ_ROOT / "bin"
ext = ".exe" if platform.system() == "Windows" else ""
exe = f"mf6{ext}"


def test_cli_version():
    output = " ".join(
        subprocess.check_output([str(bin_path / exe), "-v"]).decode().split()
    )
    print(output)
    assert output.startswith(f"{exe}:"), f"found: {output}"
    version = output.lower().split(" ")[1]
    print(version)
    v_split = version.split(".")
    assert len(v_split) >= 2
    assert all(s[-1].isdigit() for s in v_split[:2])
