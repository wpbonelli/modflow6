from itertools import chain
from os import PathLike
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from modflow_devtools.dfn import Dfn

from filters import Filters


def _get_template_env():
    template_loader = FileSystemLoader(Path(__file__).parent)
    template_env = Environment(
    loader=template_loader,
    trim_blocks=True,
    lstrip_blocks=True,
    line_statement_prefix="_",
    keep_trailing_newline=True,
    )
    template_env.filters["value"] = Filters.value
    return template_env


def make_selectors(dfns: dict[str, Dfn], outdir: PathLike, verbose: bool = False):
    pass


def make_components(dfn: Dfn, outdir: PathLike, verbose: bool = False):
    env = _get_template_env()
    outdir = Path(outdir).expanduser().resolve().absolute()
    template_name = "IdmComponent.idm.f90.jinja"
    template = env.get_template(template_name)
    components = list(
        chain.from_iterable(IdmComponentDescriptor.from_dfn(dfn) for dfn in dfns.values())
    )
    for component in components:
        component_name = component["name"]
        target_path = outdir / f"{Filters.title(component_name)}.f90"
        with open(target_path, "w") as f:
            f.write(template.render(**component))
            if verbose:
                print(f"Wrote {target_path}")


def make_all(dfndir: PathLike, outdir: PathLike, verbose: bool = False, version: int = 1):
    """Generate Fortran source files from DFN files."""

    dfndir = Path(dfndir).expanduser().resolve().absolute()
    dfns = Dfn.load_all(dfndir, version=version)
    make_selectors(dfns, outdir, verbose)
    for dfn in dfns.values():
        make_components(dfn, outdir, verbose)