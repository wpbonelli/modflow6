"""Write one fixed MF6 input set using the same config as
test_mmr_to_mf6_dfw_regression. Args: domain_dir inflow_dir run_dir"""

import pathlib as pl
import sys

import pywatershed as pws
from pywatershed.utils.mmr_to_mf6_dfw import MmrToMf6Dfw

domain_dir, inflow_dir, run_dir = map(pl.Path, sys.argv[1:4])

seg_shp_file = pws.utils.gis_files.gis_dir / "drb_2yr/Segments_subset.shp"
control = pws.Control.load_prms(
    domain_dir / "nhm.control", warn_unused_options=False
)
ndays_run = 18
control.edit_end_time(
    control.start_time + ((ndays_run - 1) * control.time_step)
)
params = pws.Parameters.merge(
    pws.Parameters.from_netcdf(domain_dir / "parameters_dis_both.nc"),
    pws.Parameters.from_netcdf(domain_dir / "parameters_PRMSChannel.nc"),
)
hclose, rclose, relax = 1e-5, 1.0, 0.97
ims_options = {
    "print_option": "SUMMARY",
    "outer_dvclose": hclose,
    "outer_maximum": 100,
    "under_relaxation": "DBD",
    "under_relaxation_theta": 0.95,
    "under_relaxation_kappa": 0.0001,
    "under_relaxation_gamma": 0.0,
    "under_relaxation_momentum": 0.0,
    "inner_maximum": 50,
    "inner_dvclose": hclose,
    "linear_acceleration": "BICGSTAB",
    "scaling_method": "NONE",
    "reordering_method": "NONE",
    "relaxation_factor": relax,
    "filename": "drb.ims",
}
run_dir.mkdir(parents=True, exist_ok=True)
dfw = MmrToMf6Dfw(
    control=control,
    segment_shp_file=seg_shp_file,
    params=params,
    tdis_perlen=2 * 3600,
    tdis_nstp=3,
    output_dir=run_dir,
    sim_name="drb_dfw",
    ims_options=ims_options,
    dfw_options={"print_flows": True, "save_flows": True, "idcxs": None},
    sto_options={"save_flows": True},
    oc_options={
        "saverecord": [("STAGE", "ALL"), ("BUDGET", "ALL")],
        "printrecord": [("STAGE", "LAST"), ("BUDGET", "ALL")],
    },
    ic_options={"strt": 0.5},
    chd_options={"print_input": True, "print_flows": False},
    bc_binary_files=True,
    bc_flows_combined=True,
    inflow_dir=inflow_dir,
)
dfw.write()
print("wrote", run_dir)

# flopy 3.10.0 drops CHD cellids; restore them (values/format untouched).
import flopy
if flopy.__version__ != "3.10.0":
    print("flopy", flopy.__version__, "no CHD patch")
    sys.exit(0)
chd = run_dir / "drb_dfw.chd"
txt = chd.read_text().splitlines()
spd = list(dfw._outlet_chds.items())
i0 = txt.index("BEGIN period  1") + 1
assert len(txt[i0:i0 + len(spd)]) == len(spd)
for ii, (node, val) in enumerate(spd):
    line = txt[i0 + ii].strip()
    assert abs(float(line) - val) < 1e-6 * abs(val), (line, val)
    txt[i0 + ii] = f"  {node + 1}  {line}"
chd.write_text("\n".join(txt) + "\n")
