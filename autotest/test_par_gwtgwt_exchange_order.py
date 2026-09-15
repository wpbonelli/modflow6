"""
Reuses the simulation data in test_gwtgwt_exchange_order and runs it in
parallel on two processes, so that gwfa/gwta and gwfb/gwtb each land on
separate processes.

same_exchange_cells only compares a side of the exchange pair where both
models are local to the current process. With the GWT-GWT exchange's model
order reversed relative to the GWF-GWF exchange (as this scenario sets up),
neither side is local on either process, so neither process can validate the
match -- unlike in the serial case, where every model is local and the
mismatched second connection is always caught. If the match is nonetheless
accepted without comparing any cells, the expected "Cannot find GWF-GWF
exchange" error will not appear.
"""

import pytest
from framework import TestFramework
from test_gwtgwt_exchange_order import build_models as build
from test_gwtgwt_exchange_order import cases as serial_cases
from test_gwtgwt_exchange_order import check_output as check

cases = ["par_" + name for name in serial_cases]


def build_models(idx, test):
    return build(idx, test)


def check_output(idx, test):
    check(idx, test)


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        parallel=True,
        ncpus=2,
        xfail=True,
    )
    test.run()
