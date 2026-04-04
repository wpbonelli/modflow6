"""
Test a PRT model exchange-coupled to an ATS-enabled GWF model.
"""

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from flopy.utils.binaryfile import CellBudgetFile, HeadFile
from framework import TestFramework
from prt_test_utils import get_model_name
from test_prt_exg import build_mf6_sim

simname = "prt_ats"
cases = [simname]

# ATS parameters
dt0 = 0.5
dtmin = 1.0e-3
dtmax = 1.0
dtadj = 2.0
dtfailadj = 5.0


def build_sim(name, ws, mf6):
    sim = build_mf6_sim(name, ws, mf6)

    # add ATS to GFW TDIS
    tdis = sim.get_package("tdis")
    ats_filerecord = f"{name}.ats"
    atsperiod = [(0, dt0, dtmin, dtmax, dtadj, dtfailadj)]
    tdis.ats.initialize(
        maxats=len(atsperiod),
        perioddata=atsperiod,
        filename=ats_filerecord,
    )

    return sim


def build_models(idx, test):
    return build_sim(test.name, test.workspace, test.targets["mf6"])


def check_output(idx, test):
    name = test.name
    ws = test.workspace
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    gwf_head_file = ws / f"{gwf_name}.hds"
    gwf_budget_file = ws / f"{gwf_name}.bud"

    gwf_hds = HeadFile(gwf_head_file)
    gwf_cbb = CellBudgetFile(gwf_budget_file)
    gwf_times = set(gwf_hds.get_times())
    assert gwf_times == {0.5, 1.0}
    assert gwf_times == set(gwf_cbb.get_times())

    prt_track_csv_file = ws / f"{prt_name}.trk.csv"
    prt_budget_file = ws / f"{prt_name}.cbb"

    prt_pls = pd.read_csv(prt_track_csv_file)
    prt_cbb = CellBudgetFile(prt_budget_file)
    assert gwf_times == set(prt_cbb.get_times())

    timestep_events = prt_pls[prt_pls.ireason == 2]
    assert gwf_times == set(timestep_events.t.unique())


def plot_output(idx, test):
    name = test.name
    ws = test.workspace
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    sim = test.sims[0]
    gwf = sim.get_model(gwf_name)
    mg = gwf.modelgrid

    # check mf6 output files exist
    gwf_head_file = ws / f"{gwf_name}.hds"
    prt_track_csv_file = ws / f"{prt_name}.trk.csv"
    prt_budget_file = ws / f"{prt_name}.cbb"

    # extract head, budget, and specific discharge results from GWF model
    hds = HeadFile(ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    prt_pls = pd.read_csv(ws / prt_track_csv_file, na_filter=False)
    prt_bud = flopy.utils.CellBudgetFile(prt_budget_file, precision="double")

    # set up plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    ax.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax)
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mf6_plines = prt_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            title="MF6 pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax,
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
        )

    # view/save plot
    plt.show()
    plt.savefig(ws / f"{name}.png")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
    )
    test.run()
