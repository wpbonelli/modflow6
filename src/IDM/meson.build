# paired Fortran files for DFNs
py = find_program('python', required : true) 
if py.found()
  # todo: once IDM fully integrated, discover all DFN files automatically?
  dfns = fs.read(root_path / 'utils' / 'idmloader' / 'dfns.txt').split('\n')
  args = [
      root_path / 'utils' / 'idmloader' / 'scripts' / 'dfn2f90.py',
      '-d', '@INPUT@', '-o', '@BUILD_DIR@',
      # '-v'  # uncomment for verbose
    ]
  gen = generator(py, output: '@BASENAME@idm.f90', arguments: args)
  foreach dfn : dfns
    if dfn.startswith('#') or dfn.strip() == ''
      continue  # skip comments and blank lines
    endif
    modflow_sources += gen.process(dfn_path / dfn)
  endforeach
else
  error('Can\'t find python, can\'t generate IDM source files')
endif

# definition selector files
modflow_sources += files(
  'IdmDfnSelector.f90',
  'IdmExgDfnSelector.f90',
  'IdmSimDfnSelector.f90')
foreach m : models
  mcc = m[0].to_upper() + m.substring(1, 3)  # camel case
  modflow_sources += files(f'Idm@mcc@DfnSelector.f90')
endforeach