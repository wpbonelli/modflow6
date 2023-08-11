import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from shapely.geometry import MultiPoint, LineString
import flopy
from flopy.utils.gridintersect import GridIntersect
import numpy as np
from pathlib import Path

sim_name = "tetratech"
nparts = 2
sim_ws = Path("temp") / f"mf6-orig"
spl_ws = Path("temp") / f"mf6-spl-{nparts}"
sim_ws.mkdir(exist_ok=True, parents=True)
spl_ws.mkdir(exist_ok=True, parents=True)

from flopy.mf6 import MFSimulation
sim = MFSimulation.load(sim_name, sim_ws=sim_ws)#, load_only=["dis", "npf", "ic"])

# Generate the splitting array.

from flopy.mf6.utils import Mf6Splitter

spl = Mf6Splitter(sim)
array = spl.optimize_splitting_mask(nparts=nparts)

# View the splitting array.

fig, ax = plt.subplots(figsize=(10, 10))
ax.set_aspect("equal")
gwf = sim.get_model("carrier")
pmv = flopy.plot.PlotMapView(model=gwf, ax=ax)
pa = pmv.plot_array(array)
plt.colorbar(pa, shrink=0.5)
plt.show()

# Split the model.

spl_sim = spl.split_model(array)
spl_sim.set_sim_path(spl_ws)
spl_sim.write_simulation()