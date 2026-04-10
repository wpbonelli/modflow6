# PRT Porosity Validation Implementation Plan

## Objective
Add porosity validation to the PRT (Particle Tracking) model to trap invalid porosity values and handle edge cases properly.

## Requirements
1. Trap negative porosity in the MIP (Model Input Package) - prevent simulation from running
2. Validate that porosity ≤ 1.0 (physically realistic upper bound)
3. Handle zero porosity - two approaches documented below (decision deferred)

## Background

### Current Issue
- Porosity is read in `prt-mip.f90` but NOT validated for negative, zero, or > 1.0 values
- Zero or negative porosity causes division by zero in velocity calculations:
  - `CellUtil.f90:100`: `factor = DONE / (defn%retfactor * defn%porosity)`
  - `MethodDis.f90:120`: `factor = factor / cell%defn%porosity`
- Currently only checks that porosity is specified, not that values are valid

### Particle Status Codes
- `TERM_NO_EXITS = 5`: Terminated in cell with no exit face (appropriate for zero porosity)
- Defined in `Particle.f90:29-39`

## Implementation Approaches

### Approach A: Reject Zero Porosity at Input Time (Simpler)

**Validation in prt-mip.f90 only:**
```fortran
! In mip_ar subroutine, after line 137
do n = 1, this%dis%nodes
  if (this%porosity(n) <= DZERO) then
    write (errmsg, '(a,i0,a,g0,a)') &
      'POROSITY for cell ', n, &
      ' must be greater than 0.0 (specified value is ', &
      this%porosity(n), ').'
    call store_error(errmsg)
  else if (this%porosity(n) > DONE) then
    write (errmsg, '(a,i0,a,g0,a)') &
      'POROSITY for cell ', n, &
      ' must be less than or equal to 1.0 (specified value is ', &
      this%porosity(n), ').'
    call store_error(errmsg)
  end if
end do

if (count_errors() > 0) then
  call store_error_filename(this%input_fname)
end if
```

**Pros:**
- Simpler implementation (single file change)
- Fails fast with clear error messages
- No runtime overhead during tracking
- Follows MODFLOW 6 "validate early" pattern

**Cons:**
- Rejects models with zero porosity cells
- Less flexible for special use cases

**Files to modify:**
- `src\Model\ParticleTracking\prt-mip.f90`

### Approach B: Allow Zero Porosity, Check at Tracking Time

**Validation in prt-mip.f90:**
```fortran
! Same as Approach A, but check: if (this%porosity(n) < DZERO)
! Allow zero, reject only negative values
```

**Additional checks before velocity calculations:**

**In CellUtil.f90** (before line 100):
```fortran
! Check for zero porosity to avoid division by zero
if (defn%porosity <= DZERO) then
  rect%vx1 = DZERO
  rect%vx2 = DZERO
  rect%vy1 = DZERO
  rect%vy2 = DZERO
  rect%vz1 = DZERO
  rect%vz2 = DZERO
  defn%inoexitface = 1  ! Flag as no exit face
  return
end if
```

**In MethodDis.f90** (before line 119):
```fortran
! Check for zero porosity
if (cell%defn%porosity <= DZERO) then
  cell%vx1 = DZERO
  cell%vx2 = DZERO
  cell%vy1 = DZERO
  cell%vy2 = DZERO
  cell%vz1 = DZERO
  cell%vz2 = DZERO
  cell%defn%inoexitface = 1
  return
end if
```

**In MethodDisv.f90** (similar location in load_cell):
```fortran
! Same check as MethodDis.f90
```

**Pros:**
- Allows zero porosity in input files
- Particles terminate gracefully with TERM_NO_EXITS status
- More flexible for edge cases

**Cons:**
- More complex (3-4 file changes)
- Small runtime overhead for porosity checks during tracking
- Division by zero protection must be added in multiple locations

**Files to modify:**
- `src\Model\ParticleTracking\prt-mip.f90`
- `src\Solution\ParticleTracker\Domain\CellUtil.f90`
- `src\Solution\ParticleTracker\Method\MethodDis.f90`
- `src\Solution\ParticleTracker\Method\MethodDisv.f90`

## Common Changes (Both Approaches)

### 1. Update imports in prt-mip.f90
Line 11, change:
```fortran
use SimModule, only: store_error
```
to:
```fortran
use SimModule, only: store_error, store_error_filename, count_errors
```

### 2. Add validation loop in mip_ar subroutine
After line 137 in `prt-mip.f90`, add porosity validation loop as shown in each approach.

## Critical Files

1. **`src\Model\ParticleTracking\prt-mip.f90`** (Required for both approaches)
   - Line 11: Update imports
   - Lines 133-138: Modify mip_ar subroutine

2. **`src\Solution\ParticleTracker\Domain\CellUtil.f90`** (Approach B only)
   - Line 100: Add zero porosity check before division

3. **`src\Solution\ParticleTracker\Method\MethodDis.f90`** (Approach B only)
   - Line 119-120: Add zero porosity check before division

4. **`src\Solution\ParticleTracker\Method\MethodDisv.f90`** (Approach B only)
   - Similar location: Add zero porosity check

## Validation Pattern

Following MODFLOW 6 conventions from `gwf-csub.f90`, `gwf-uzf.f90`:
- Loop through all nodes/cells
- Check each porosity value
- Use `write(errmsg, ...)` + `call store_error(errmsg)` to accumulate errors
- After loop: `if (count_errors() > 0) call store_error_filename(this%input_fname)`
- This reports ALL invalid cells in one run, not just the first error

## Testing Recommendations

1. **Negative porosity**: Should fail with clear error message
2. **Zero porosity**:
   - Approach A: Should fail at input time
   - Approach B: Should run, particles terminate with TERM_NO_EXITS
3. **Porosity > 1.0**: Should fail with clear error message
4. **Valid range (0 < porosity ≤ 1.0)**: Should run successfully
5. **Multiple invalid cells**: Should report all errors together

## Decision Point

The implementation will follow **Approach A or B** based on final decision about zero porosity handling. Both approaches validate that porosity ≤ 1.0 and reject negative values.
