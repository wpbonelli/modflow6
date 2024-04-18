from pathlib import Path
from tempfile import TemporaryDirectory

from autotest.conftest import get_targets
from autotest.framework import TestFramework


class PrtSuite:

    def time_prt_triangle(self):
        from autotest.test_prt_triangle import build_models, check_output

        with TemporaryDirectory() as tmpdir:
            test = TestFramework(
                name="prt-triangle",
                workspace=tmpdir.name,
                build=lambda t: build_models(0, t),
                check=lambda t: check_output(0, t),
                targets=get_targets(),
                compare=None,
            )
            test.run()

    def time_prt_voronoi1(self):
        from autotest.test_prt_voronoi1 import build_models, check_output

        with TemporaryDirectory() as tmpdir:
            test = TestFramework(
                name="prt-voronoi1",
                workspace=tmpdir.name,
                build=lambda t: build_models(0, t),
                check=lambda t: check_output(0, t),
                targets=get_targets(),
                compare=None,
            )
            test.run()
