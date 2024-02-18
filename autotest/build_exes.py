import argparse
from pathlib import Path

import pytest
from modflow_devtools.build import meson_build

from conftest import PROJ_ROOT

repository = "MODFLOW-USGS/modflow6"
top_bin_path = PROJ_ROOT / "bin"


@pytest.fixture
def bin_path():
    return top_bin_path


def test_meson_build(bin_path):
    meson_build(
        project_path=PROJ_ROOT,
        build_path=PROJ_ROOT / "builddir",
        bin_path=bin_path,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Rebuild local development version of MODFLOW 6"
    )
    parser.add_argument(
        "-p", "--path", help="path to bin directory", default=top_bin_path
    )
    args = parser.parse_args()
    test_meson_build(Path(args.path).expanduser().resolve())
