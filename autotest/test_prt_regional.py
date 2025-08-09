from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import flopy


PROJ_ROOT = Path(__file__).parents[1]


def test_prt_regional(targets, plot):
    ws = PROJ_ROOT / "autotest" / "temp" / "regional_PRT_pre_ag_IFLOWFACE"
    sim = flopy.mf6.MFSimulation.load(
        sim_ws=ws,
        version="mf6",
        exe_name=targets["mf6"],
    )
    # sim.run_simulation(silent=False)
    gwf_name = "regional"
    prt_name = "reg_prt"
    gwf = sim.get_model(gwf_name)
    mg = gwf.modelgrid
    prt_track_csv_file = ws / f"{prt_name}.csv"
    mf6_pls = pd.read_csv(prt_track_csv_file, na_filter=False)
    endpt = mf6_pls[mf6_pls.istatus > 1].iloc[-1]
    nodes = [nn - 1 for nn in [201564, 201563, 201030, 505906, 506439]]
    xc, yc = mg.get_xcellcenters_for_layer(0), mg.get_ycellcenters_for_layer(0)
    cellids = mg.get_lrc(nodes)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    ax[0].set_aspect("equal")
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax[0], layer=0)
    pmv.plot_grid(lw=0.1)
    pmv.plot_bc("SFR", plotAll=True, alpha=0.2)
    pmv.plot_bc("GHB", plotAll=True, alpha=0.2)
    row = None
    column = None
    for nn, (k, i, j) in zip(nodes, cellids):
        x, y = xc[i, j], yc[i, j]
        if row is None:
            row = i
            column = j - 1
        ax[0].plot(x, y, "o", color="grey", alpha=0.5)
        ax[0].annotate(f"{i + 1}, {j + 1}", (x, y), color="grey", alpha=0.5)
        vertices = mg.get_cell_vertices(node=nn)
        vertices_with_z = [(vx, vy, mg.zvertices[k, i, j]) for vx, vy in vertices]
        for vx, vy, vz in vertices_with_z:
            ax[0].plot(vx, vy, "o", color="blue", alpha=0.5)
    mf6_pls.plot(x="x", y="y", ax=ax[0], color="black", marker=".", markersize=3, alpha=0.5, label=None)
    
    hpad = 1000
    ax[0].set_xlim(left=endpt.x - hpad, right=endpt.x + hpad)
    ax[0].set_ylim(bottom=endpt.y - hpad, top=endpt.y + hpad)
    ax[0].set_title("Map view")

    ax[1].set_aspect(10)
    pxs1 = flopy.plot.PlotCrossSection(model=gwf, ax=ax[1], line={"row": row})
    pxs1.plot_grid(lw=0.1)
    pxs1.plot_bc("SFR", alpha=0.2)
    pxs1.plot_bc("GHB", alpha=0.2)
    for nn, (k, i, j) in zip(nodes, cellids):
        if i != row:
            continue
        x, y = xc[i, j], mg.zcellcenters[k, i, j]
        ax[1].plot(x, y, "o", color="grey", alpha=0.5)
        ax[1].annotate(f"{k + 1}, {i + 1}, {j + 1}", (x, y), color="grey")
        vertices = mg.get_cell_vertices(node=nn)
        vertices_with_z = [(vx, vy, mg.zvertices[k, i, j]) for vx, vy in vertices]
        for vx, vy, vz in vertices_with_z:
            ax[1].plot(vx, vz, "o", color="blue", alpha=0.5)
    vpad = 50
    ax[1].set_xlim(left=endpt.x - hpad, right=endpt.x + hpad)
    ax[1].set_ylim(bottom=endpt.z - vpad, top=endpt.z + vpad)
    ax[1].set_title(f"Row {row + 1}")
    mf6_pls.plot(x="x", y="z", ax=ax[1], color="black", marker=".", markersize=3, alpha=0.5, label=None)

    # TODO why grid wrong?
    # ax[2].set_aspect(10)
    # pxs2 = flopy.plot.PlotCrossSection(model=gwf, ax=ax[2], line={"column": column})
    # pxs2.plot_grid(lw=0.1)
    # pxs2.plot_bc("SFR", alpha=0.2)
    # pxs2.plot_bc("GHB", alpha=0.2)
    # pxs2.plot_pathline(mf6_pls, method="cell")
    # for nn, (k, i, j) in zip(nodes, cellids):
    #     if j != column:
    #         continue
    #     x, y = yc[i, j], mg.zcellcenters[k, i, j]
    #     ax[2].plot(x, y, "o", color="grey", alpha=0.5)
    #     ax[2].annotate(f"{k + 1}, {i + 1}, {j + 1}", (x, y), color="grey", alpha=0.5)
    #     vertices = mg.get_cell_vertices(node=nn)
    #     vertices_with_z = [(vx, vy, mg.zvertices[k, i, j]) for vx, vy in vertices]
    #     for vx, vy, vz in vertices_with_z:
    #         ax[2].plot(vy, vz, "o", color="blue", alpha=0.5)
    # ax[2].set_xlim(left=endpt.y - hpad, right=endpt.y + hpad)
    # ax[2].set_ylim(bottom=endpt.z - vpad, top=endpt.z + vpad)
    # ax[2].set_title(f"Column {column + 1}")
    # mf6_pls.plot(x="y", y="z", ax=ax[2], color="black", marker=".", markersize=3, alpha=0.5, label=None)

    plt.show()
