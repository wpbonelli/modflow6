"""Build synthetic models that stress each TDIS drift mechanism."""

import sys
from pathlib import Path

import flopy

out = Path(sys.argv[1])

CASES = {
    # name: (perioddata, time_units, ats)
    "uniform_36525": ([(36525.0, 36525, 1.0)], "days", False),
    "uniform_1e5_tenth": ([(1.0e4, 100000, 1.0)], "days", False),
    "geom_r1.01_n1000": ([(1000.0, 1000, 1.01)], "days", False),
    "geom_r1.05_n500": ([(1.0e5, 500, 1.05)], "days", False),
    "geom_r1.2_n50": ([(100.0, 50, 1.2)], "days", False),
    "shrink_r0.95_n100": ([(100.0, 100, 0.95)], "days", False),
    "shrink_r0.99_n1000": ([(100.0, 1000, 0.99)], "days", False),
    "shrink_r0.5_n50": ([(100.0, 50, 0.5)], "days", False),
    "near1_r1.0001_n1000": ([(100.0, 1000, 1.0001)], "days", False),
    "near1_r1.00001_n10000": ([(100.0, 10000, 1.00001)], "days", False),
    "near1_r0.9999_n1000": ([(100.0, 1000, 0.9999)], "days", False),
    "near1_r1.003_n1000": ([(100.0, 1000, 1.003)], "days", False),
    "periods_2000x0.1": ([(0.1, 1, 1.0)] * 2000, "days", False),
    "periods_1200x_month_r1.1": ([(30.4375, 10, 1.1)] * 1200, "days", False),
    "periods_1000x_r0.9": ([(1.0, 10, 0.9)] * 1000, "days", False),
    "seconds_100yr_daily": ([(3.15576e9, 36525, 1.0)], "seconds", False),
    "seconds_100x_year_r1.05": ([(3.15576e7, 365, 1.05)] * 100, "seconds", False),
    # ATS: (dt0, dtmax); retries are forced with MF6_TDIS_FORCE_RETRY
    "ats_50x10.1": ([(10.1, 1, 1.0)] * 50, "days", (1.0, 3.0)),
    "ats_seconds_100x_year": ([(3.15576e7, 1, 1.0)] * 100, "seconds", (1.0e6, 5.0e6)),
}


def build(name, pd, units, ats):
    ws = out / name
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, verbosity_level=0)
    tdis = flopy.mf6.ModflowTdis(sim, time_units=units, nper=len(pd), perioddata=pd)
    if ats:
        dt0, dtmax = ats
        ap = [(i, dt0, 1.0e-5, dtmax, 2.0, 5.0) for i in range(len(pd))]
        tdis.ats.initialize(maxats=len(ap), perioddata=ap, filename=f"{name}.ats")
    nlay, botm = 1, [0.0]
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    ims = {"ats_outer_maximum_fraction": 0.5} if ats else {}
    flopy.mf6.ModflowIms(sim, outer_maximum=5, inner_maximum=10, **ims)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=1, ncol=1, delr=1.0, delc=1.0, top=100.0, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=90.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1e-5, transient={0: True})
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data={0: [[(0, 0, 0), 90.0, 1.0]]})
    sim.write_simulation(silent=True)


names = sys.argv[2:] or list(CASES)
for n in names:
    build(n, *CASES[n])
    print("built", n)
