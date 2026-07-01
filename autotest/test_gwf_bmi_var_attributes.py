"""
Test that the memory manager's readonly/output variable attributes are
enforced uniformly through the BMI, not just by convention:

- get_input_var_names()/get_output_var_names() are filtered on these
  attributes rather than dumping the whole memory store.
- a readonly variable (WEL's SIMVALS) cannot be written via set_value(),
  regardless of which BMI entry point is used.
- a variable that is both writable and output (the model's X/head array)
  can legitimately appear in both lists and be written.
"""

import flopy
import numpy as np
import pytest
from modflow_devtools.markers import requires_pkg

name = "bmiattrs"
nlay, nrow, ncol = 1, 1, 5
delr = delc = 1.0
top = 1.0
botm = [0.0]
k = 1.0
strt = 10.0
chd_head = 10.0
wellq = -1.0


@pytest.fixture
def simple_sim(tmp_path):
    """A trivial, linear, steady-state 1D column: CHD on one end, WEL on the other."""
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(tmp_path))
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-8, inner_dvclose=1e-8)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), chd_head]])
    flopy.mf6.ModflowGwfwel(
        gwf, pname="WEL_0", stress_period_data=[[(0, 0, ncol - 1), wellq]]
    )
    sim.write_simulation()
    return sim


@requires_pkg("xmipy")
def test_var_attributes(simple_sim, targets):
    from xmipy import XmiWrapper
    from xmipy.errors import XMIError

    sim = simple_sim
    mf6 = XmiWrapper(lib_path=targets["libmf6"], working_directory=sim.sim_path)
    mf6.initialize()

    simvals_addr = mf6.get_var_address("SIMVALS", name, "WEL_0")
    x_addr = mf6.get_var_address("X", name)

    input_vars = mf6.get_input_var_names()
    output_vars = mf6.get_output_var_names()

    # SIMVALS is readonly + output: input-excluded, output-included
    assert simvals_addr not in input_vars
    assert simvals_addr in output_vars

    # X is writable + output: appears in both
    assert x_addr in input_vars
    assert x_addr in output_vars

    # writing SIMVALS should be rejected uniformly, regardless of entry point
    simvals = mf6.get_value_ptr(simvals_addr)
    with pytest.raises(XMIError, match="read-only"):
        mf6.set_value(simvals_addr, np.zeros_like(simvals))

    # writing X should succeed, since it's not readonly
    x = mf6.get_value_ptr(x_addr)
    new_x = np.full_like(x, 5.0)
    mf6.set_value(x_addr, new_x)
    np.testing.assert_array_equal(mf6.get_value_ptr(x_addr), new_x)

    mf6.finalize()
