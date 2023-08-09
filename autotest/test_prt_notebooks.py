import re
from os import environ
from pathlib import Path
from platform import system
from pprint import pprint

import pytest
from flaky import flaky
from modflow_devtools.markers import excludes_platform
from modflow_devtools.misc import run_cmd, set_env

from conftest import project_root_path


def get_notebook_scripts(pattern=None, exclude=None):
    repos_path = environ.get("REPOS_PATH", None)
    if repos_path is None:
        repos_path = project_root_path.parent
    repo_path = Path(repos_path) / "modflow6-examples"
    if not repo_path.is_dir():
        return []
    nbpaths = [
        str(p)
        for p in (repo_path / "scripts").glob("*.py")
        if pattern is None or pattern in p.name
    ]

    # sort for pytest-xdist: workers must collect tests in the same order
    return sorted(
        [p for p in nbpaths if not exclude or not any(e in p for e in exclude)]
    )


# @flaky(max_runs=3)
# @excludes_platform("Windows", ci_only=True)
@pytest.mark.slow
@pytest.mark.parametrize(
    "notebook",
    get_notebook_scripts(pattern="ex-prt", exclude=["ex-prt-mp7-p03"]),
)
def test_notebooks(notebook, function_tmpdir, targets):
    delim = ";" if system() == "Windows" else ":"
    path = (
        environ.get("PATH", "")
        + f"{delim}{targets.mf6.parent}"
        + f"{delim}{targets.mf6.parent / 'downloaded'}"
        + f"{delim}{targets.mf6.parent / 'rebuilt'}"
    )
    with set_env(PATH=path):
        args = [
            "jupytext",
            "--from",
            "py",
            "--to",
            "ipynb",
            "--execute",
            "--run-path",
            function_tmpdir,
            "--output",
            function_tmpdir / f"{Path(notebook).stem}.ipynb",
            notebook,
        ]
        stdout, stderr, returncode = run_cmd(*args, verbose=True)

    if returncode != 0:
        if "Missing optional dependency" in stderr:
            pkg = re.findall("Missing optional dependency '(.*)'", stderr)[0]
            pytest.skip(f"notebook requires optional dependency {pkg!r}")
        elif "No module named " in stderr:
            pkg = re.findall("No module named '(.*)'", stderr)[0]
            pytest.skip(f"notebook requires package {pkg!r}")

    assert returncode == 0, f"could not run {notebook}"
    pprint(stdout)
    pprint(stderr)
