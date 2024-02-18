import sys
from pathlib import Path
from typing import Dict
from warnings import warn

import pytest
from syrupy.extensions.single_file import SingleFileSnapshotExtension
from modflow_devtools.ostags import get_binary_suffixes

pytest_plugins = ["modflow_devtools.fixtures"]  # needs to be lowercase
PROJ_ROOT = Path(__file__).resolve().parent.parent
EXE_EXT, LIB_EXT = get_binary_suffixes(sys.platform)
BIN_PATH = PROJ_ROOT / "bin"
BIN_PATH_DL = BIN_PATH / "downloaded"
BIN_PATH_RB = BIN_PATH / "rebuilt"
BINARIES = {
    "development": [
        ("mf6", BIN_PATH / f"mf6{EXE_EXT}"),
        ("libmf6", BIN_PATH / f"libmf6{LIB_EXT}"),
        ("mf5to6", BIN_PATH / f"mf5to6{EXE_EXT}"),
        ("zbud6", BIN_PATH / f"zbud6{EXE_EXT}"),
    ],
    "downloaded": [
        ("mf2000", BIN_PATH_DL / f"mf2000{EXE_EXT}"),
        ("mf2005", BIN_PATH_DL / f"mf2005dbl{EXE_EXT}"),
        ("mfnwt", BIN_PATH_DL / f"mfnwtdbl{EXE_EXT}"),
        ("mfusg", BIN_PATH_DL / f"mfusgdbl{EXE_EXT}"),
        ("mflgr", BIN_PATH_DL / f"mflgrdbl{EXE_EXT}"),
        ("mf2005s", BIN_PATH_DL / f"mf2005{EXE_EXT}"),
        ("mt3dms", BIN_PATH_DL / f"mt3dms{EXE_EXT}"),
        ("crt", BIN_PATH_DL / f"crt{EXE_EXT}"),
        ("gridgen", BIN_PATH_DL / f"gridgen{EXE_EXT}"),
        ("mp6", BIN_PATH_DL / f"mp6{EXE_EXT}"),
        ("mp7", BIN_PATH_DL / f"mp7{EXE_EXT}"),
        ("swtv4", BIN_PATH_DL / f"swtv4{EXE_EXT}"),
        ("sutra", BIN_PATH_DL / f"sutra{EXE_EXT}"),
        ("triangle", BIN_PATH_DL / f"triangle{EXE_EXT}"),
        ("vs2dt", BIN_PATH_DL / f"vs2dt{EXE_EXT}"),
        ("zonbudusg", BIN_PATH_DL / f"zonbudusg{EXE_EXT}"),
    ],
    "rebuilt": [
        ("mf6_regression", BIN_PATH_RB / f"mf6{EXE_EXT}"),
        ("libmf6_regression", BIN_PATH_RB / f"libmf6{LIB_EXT}"),
        ("mf5to6_regression", BIN_PATH_RB / f"mf5to6{EXE_EXT}"),
        ("zbud6_regression", BIN_PATH_RB / f"zbud6{EXE_EXT}"),
    ],
}


@pytest.fixture(scope="session")
def bin_path() -> Path:
    return BIN_PATH


@pytest.fixture(scope="session")
def targets() -> Dict[str, Path]:
    """
    Target executables for tests. These include local development builds as
    well as binaries 1) downloaded from GitHub and 2) rebuilt from the last
    official release.
    """

    d = dict()
    for k, v in BINARIES["development"]:
        # require development binaries
        assert v.is_file(), f"Couldn't find binary '{k}' expected at: {v}"
        d[k] = v
    for k, v in BINARIES["downloaded"]:
        # downloaded binaries are optional
        if v.is_file():
            d[k] = v
        else:
            warn(f"Couldn't find downloaded binary '{k}' expected at: {v}")
    for k, v in BINARIES["rebuilt"]:
        # rebuilt binaries are optional
        if v.is_file():
            d[k] = v
        else:
            warn(f"Couldn't find rebuilt binary '{k}' expected at: {v}")
    return d


def try_get_target(targets: Dict[str, Path], name: str) -> Path:
    """Try to retrieve the path to a binary. If the binary is a development
    target and can't be found, an error is raised. Otherwise (if the binary
    is downloaded or rebuilt) the test is skipped. This is to allow testing
    without downloaded or rebuilt binaries, e.g. if the network is down."""

    exe = targets.get(name)
    if exe:
        return exe
    elif name in BINARIES["development"]:
        raise ValueError(f"Couldn't find binary '{name}'")
    else:
        pytest.skip(f"Couldn't find binary '{name}'")


@pytest.fixture
def original_regression(request) -> bool:
    return request.config.getoption("--original-regression")


@pytest.fixture(scope="session")
def markers(pytestconfig) -> str:
    return pytestconfig.getoption("-m")


def pytest_addoption(parser):
    parser.addoption(
        "--original-regression",
        action="store_true",
        default=False,
        help="use non-MF6 models for regression tests",
    )
    parser.addoption(
        "--parallel",
        action="store_true",
        default=False,
        help="include parallel test cases",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--parallel"):
        # --parallel given in cli: do not skip parallel tests
        return
    skip_parallel = pytest.mark.skip(reason="need --parallel option to run")
    for item in items:
        if "parallel" in item.keywords:
            item.add_marker(skip_parallel)


class CbcFileExtension(SingleFileSnapshotExtension):
    _file_extension = "cbc"


class HedFileExtension(SingleFileSnapshotExtension):
    _file_extension = "hed"


class HdsFileExtension(SingleFileSnapshotExtension):
    _file_extension = "hds"


@pytest.fixture
def cell_budget_file_snapshot(snapshot):
    return snapshot.use_extension(CbcFileExtension)


@pytest.fixture
def head_file_snapshot(snapshot):
    return snapshot.use_extension(HedFileExtension)
