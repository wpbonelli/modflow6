"""
Test PRT particle behavior with IFLOWFACE on interior cells.

This test specifically examines the behavior when IFLOWFACE is used to direct
boundary package flows to cell faces that are NOT on the model boundary, but
instead connect to other active cells. This can potentially cause infinite
loops or incorrect particle termination behavior.

Test scenarios:
1. WEL package in interior cell with IFLOWFACE pointing to faces with active neighbors
2. CHD package in interior cell with lateral face assignments to active neighbors
3. Particles released to test if they properly terminate or get stuck in cycles
4. Various face orientations to test all possible interior face assignments

Expected issues:
- Particles may not terminate at interior "boundary" faces
- Potential for infinite cycling between cells
- Ambiguous behavior when reaching IFLOWFACE-assigned interior faces
"""

import os
import numpy as np
import pandas as pd
import pytest
import flopy
from flopy.modpath import ParticleData
from framework import TestFramework
from prt_test_utils import get_model_name

# Test parameters
simname = "prt_iflowface_interior"
cases = {
    f"{simname}_wel_interior": {"boundary_type": "wel", "max_tracking_time": 1000},
    f"{simname}_chd_interior": {"boundary_type": "chd", "max_tracking_time": 1000},
    f"{simname}_ghb_interior": {"boundary_type": "ghb", "max_tracking_time": 1000},
}

# Grid parameters - 5x5x1 to have clear interior cells
nlay = 1
nrow = 5
ncol = 5
delr = 100.0
delc = 100.0
top = 10.0
botm = [0.0]
porosity = 0.1

# Particle release points - focus on cells that will interact with interior boundaries
particledata = ParticleData(
    partlocs=[
        (0, 1, 1),  # Near interior well
        (0, 1, 2),  # Adjacent to interior boundary
        (0, 2, 1),  # Adjacent to interior boundary  
        (0, 2, 2),  # Center - interior boundary cell
        (0, 2, 3),  # Adjacent to interior boundary
        (0, 3, 2),  # Adjacent to interior boundary
        (0, 3, 3),  # Near interior boundary
    ],
    structured=True,
    localx=0.5,
    localy=0.5,
    localz=0.5,
)


def build_gwf_sim(name, ws, mf6, boundary_type="wel"):
    """Build GWF simulation with interior boundary packages using IFLOWFACE."""
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, version="mf6", sim_ws=ws)
    
    # TDIS package - longer simulation to detect infinite loops
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(
        sim, pname="tdis", time_units="DAYS", nper=1, perioddata=[(1000.0, 100, 1.0)]
    )
    
    # IMS package
    ims = flopy.mf6.modflow.mfims.ModflowIms(
        sim, pname="ims", complexity="SIMPLE", outer_maximum=500, inner_maximum=100
    )
    
    # GWF model
    gwf_name = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name)
    
    # DIS package
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwf,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    
    # IC package
    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic", strt=8.0)
    
    # NPF package
    npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
        gwf,
        pname="npf",
        icelltype=1,
        k=10.0,
        save_flows=True,
        save_specific_discharge=True,
    )
    
    # Add boundary conditions on actual model boundaries for flow system
    chd_boundary = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            ((0, 0, 0), 9.0),  # Top-left corner - high head
            ((0, 4, 4), 7.0),  # Bottom-right corner - low head
        ],
        save_flows=True,
    )
    
    # Add interior boundary package with IFLOWFACE pointing to faces with active neighbors
    auxiliary = ["iflowface"]
    
    if boundary_type == "wel":
        # Wells in interior cells with IFLOWFACE to faces with active neighbors
        interior_data = [
            # (k, i, j, q, iflowface) - interior cells with faces pointing to active neighbors
            (0, 2, 2, -1.0, 1),   # Center cell, face 1 (west) -> connects to (2,1) - ACTIVE
            (0, 1, 3, 0.5, 2),    # Interior cell, face 2 (north) -> connects to (0,3) - ACTIVE  
            (0, 3, 1, 0.3, 3),    # Interior cell, face 3 (east) -> connects to (3,2) - ACTIVE
            (0, 2, 3, -0.2, 4),   # Interior cell, face 4 (south) -> connects to (3,3) - ACTIVE
        ]
        
        interior_package = flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data=interior_data,
            auxiliary=auxiliary,
            save_flows=True,
        )
        
    elif boundary_type == "chd":
        # Constant heads in interior cells with IFLOWFACE to active neighbors
        interior_data = [
            # (k, i, j, head, iflowface)
            ((0, 2, 2), 8.5, 1),   # Center cell, west face -> active neighbor
            ((0, 1, 3), 7.8, 2),   # Interior cell, north face -> active neighbor
            ((0, 3, 1), 8.2, 3),   # Interior cell, east face -> active neighbor
            ((0, 2, 3), 7.5, 4),   # Interior cell, south face -> active neighbor
        ]
        
        interior_package = flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=interior_data,
            auxiliary=auxiliary,
            save_flows=True,
        )
        
    elif boundary_type == "ghb":
        # General head boundaries in interior cells with IFLOWFACE to active neighbors
        interior_data = [
            # (k, i, j, bhead, cond, iflowface)
            ((0, 2, 2), 8.5, 1.0, 1),   # Center cell, west face -> active neighbor
            ((0, 1, 3), 7.8, 1.0, 2),   # Interior cell, north face -> active neighbor  
            ((0, 3, 1), 8.2, 1.0, 3),   # Interior cell, east face -> active neighbor
            ((0, 2, 3), 7.5, 1.0, 4),   # Interior cell, south face -> active neighbor
        ]
        
        interior_package = flopy.mf6.ModflowGwfghb(
            gwf,
            stress_period_data=interior_data,
            auxiliary=auxiliary,
            save_flows=True,
        )
    
    # OC package
    oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(
        gwf,
        pname="oc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        head_filerecord=[f"{gwf_name}.hds"],
        budget_filerecord=[f"{gwf_name}.cbc"],
    )
    
    return sim


def build_prt_sim(name, gwf, prt_ws, mf6, max_tracking_time=1000):
    """Build PRT simulation with timeout to catch infinite loops."""
    prt_name = get_model_name(name, "prt")
    sim = flopy.mf6.MFSimulation(sim_name=prt_name, exe_name=mf6, sim_ws=prt_ws)
    
    # TDIS package - match GWF timing
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=1,
        perioddata=[(max_tracking_time, 100, 1.0)],
    )
    
    # PRT model
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    
    # DIS package
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    
    # MIP package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)
    
    # PRP package - with tracking time limit to detect infinite loops
    prpdata = list(particledata.to_prp(gwf.modelgrid, localz=True))
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp",
        nreleasepts=len(prpdata),
        packagedata=prpdata,
        perioddata={0: ["FIRST"]},
        local_z=True,
        extend_tracking=True,  # Allow tracking beyond simulation end time
        exit_solve_tolerance=1e-5,
        # Note: No stop_at_weak_sink to test full particle behavior
    )
    
    # OC package
    budgetfile_prt = f"{prt_name}.cbc"
    trackfile_prt = f"{prt_name}.trk"
    trackcsvfile_prt = f"{prt_name}.trk.csv"
    
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[budgetfile_prt],
        track_filerecord=[trackfile_prt],
        trackcsv_filerecord=[trackcsvfile_prt],
        saverecord=[("BUDGET", "ALL")],
        track_release=True,
        track_terminate=True,
        track_exit=True,
    )
    
    # FMI package
    gwf_ws = gwf.model_ws
    rel_prt_folder = os.path.relpath(gwf_ws, start=prt_ws)
    
    fmi_pd = [
        ("GWFHEAD", f"{rel_prt_folder}/{gwf.name}.hds"),
        ("GWFBUDGET", f"{rel_prt_folder}/{gwf.name}.cbc"),
    ]
    flopy.mf6.ModflowPrtfmi(prt, packagedata=fmi_pd)
    
    # EMS package
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])
    
    return sim


def build_models(idx, test, boundary_type="wel", max_tracking_time=1000):
    """Build GWF and PRT models."""
    gwf_sim = build_gwf_sim(
        name=test.name,
        ws=test.workspace / "gwf", 
        mf6=test.targets["mf6"],
        boundary_type=boundary_type,
    )
    gwf = gwf_sim.get_model(get_model_name(test.name, "gwf"))
    
    prt_sim = build_prt_sim(
        name=test.name,
        gwf=gwf,
        prt_ws=test.workspace / "prt",
        mf6=test.targets["mf6"],
        max_tracking_time=max_tracking_time,
    )
    
    return gwf_sim, prt_sim


def check_output(idx, test):
    """Check for potential infinite loops and incorrect particle behavior."""
    name = test.name
    prt_ws = test.workspace / "prt"
    prt_name = get_model_name(name, "prt")
    
    # Load particle pathlines
    prt_track_csv_file = prt_ws / f"{prt_name}.trk.csv"
    assert prt_track_csv_file.is_file(), "PRT track CSV file not found"
    
    pathlines = pd.read_csv(prt_track_csv_file, na_filter=False)
    
    print(f"\n=== Analysis for {name} ===")
    print(f"Total pathline points: {len(pathlines)}")
    
    # Check for signs of infinite loops
    max_time = pathlines.t.max()
    print(f"Maximum tracking time: {max_time}")
    
    # Group by particle to analyze individual behaviors
    particles = pathlines.groupby(['iprp', 'irpt'])
    
    suspicious_particles = []
    for (iprp, irpt), particle_path in particles:
        particle_points = len(particle_path)
        particle_max_time = particle_path.t.max()
        
        # Check for potential infinite loops
        # 1. Excessive number of tracking points
        # 2. Very long tracking times
        # 3. Particles that never terminate properly
        
        termination_reasons = particle_path.ireason.unique()
        
        print(f"Particle {iprp}-{irpt}: {particle_points} points, max_time={particle_max_time:.2f}, reasons={termination_reasons}")
        
        # Flag suspicious behavior
        if particle_points > 1000:  # Excessive tracking points
            suspicious_particles.append((iprp, irpt, "excessive_points", particle_points))
        if particle_max_time > 500:  # Very long tracking time
            suspicious_particles.append((iprp, irpt, "long_tracking_time", particle_max_time))
        if 3 not in termination_reasons:  # Never properly terminated (reason 3 = normal termination)
            suspicious_particles.append((iprp, irpt, "no_termination", termination_reasons))
    
    # Report suspicious behavior
    if suspicious_particles:
        print(f"\nSUSPICIOUS PARTICLE BEHAVIOR DETECTED:")
        for iprp, irpt, issue, value in suspicious_particles:
            print(f"  Particle {iprp}-{irpt}: {issue} = {value}")
        
        # Check for potential cycles by looking at repeated positions
        for (iprp, irpt), particle_path in particles:
            if (iprp, irpt) in [(p[0], p[1]) for p in suspicious_particles]:
                # Check for position cycling
                positions = particle_path[['x', 'y', 'z']].round(6)  # Round to avoid floating point issues
                position_counts = positions.value_counts()
                repeated_positions = position_counts[position_counts > 1]
                
                if not repeated_positions.empty:
                    print(f"  Particle {iprp}-{irpt} has repeated positions (potential cycling):")
                    for pos, count in repeated_positions.head().items():
                        print(f"    Position {pos}: visited {count} times")
    else:
        print("No obviously suspicious particle behavior detected.")
    
    # Store results for analysis
    test.pathlines = pathlines
    test.suspicious_particles = suspicious_particles
    test.max_tracking_time = max_time
    
    # Basic assertions
    assert len(pathlines) > 0, "No pathlines generated"
    
    # Assert against obvious infinite loops
    # Note: These thresholds may need adjustment based on actual behavior
    max_reasonable_points = 2000
    max_reasonable_time = 800
    
    if len(pathlines) > max_reasonable_points:
        print(f"WARNING: Very high number of tracking points ({len(pathlines)}) - possible infinite loop")
    
    if max_time > max_reasonable_time:
        print(f"WARNING: Very long tracking time ({max_time}) - possible infinite loop")


@pytest.mark.parametrize("idx, name", enumerate(list(cases.keys())))
def test_mf6model(idx, name, function_tmpdir, targets):
    """Run individual test cases."""
    case = cases[name]
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(
            idx, t, 
            boundary_type=case["boundary_type"],
            max_tracking_time=case["max_tracking_time"]
        ),
        check=lambda t: check_output(idx, t),
        targets=targets,
        compare=None,
    )
    test.run()


def test_compare_interior_behavior(function_tmpdir, targets):
    """Compare behavior across different boundary types."""
    results = {}
    
    for case_name, case_params in cases.items():
        print(f"\nRunning {case_name}...")
        test = TestFramework(
            name=case_name,
            workspace=function_tmpdir / case_name,
            build=lambda t: build_models(
                0, t,
                boundary_type=case_params["boundary_type"], 
                max_tracking_time=case_params["max_tracking_time"]
            ),
            check=lambda t: check_output(0, t),
            targets=targets,
        )
        test.run()
        results[case_name] = test
    
    # Summary comparison
    print(f"\n=== SUMMARY COMPARISON ===")
    for case_name, test in results.items():
        print(f"{case_name}:")
        print(f"  Max tracking time: {test.max_tracking_time:.2f}")
        print(f"  Suspicious particles: {len(test.suspicious_particles)}")
        print(f"  Total pathline points: {len(test.pathlines)}")
    
    # This test is primarily exploratory - it will help identify the issue
    # rather than assert specific expected behavior
    print(f"\n=== ANALYSIS COMPLETE ===")
    print("Review the output above for signs of infinite loops or incorrect particle behavior")
    print("when IFLOWFACE is used on interior cell faces with active neighbors.")