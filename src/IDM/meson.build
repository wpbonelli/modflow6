# IDM framework files (not auto-generated)
modflow_sources += files(
    'BoundInputContext.f90',
    'DefinitionSelect.f90',
    'IdmLoad.f90',
    'IdmLogger.f90',
    'InputDefinition.f90',
    'InputLoadType.f90',
    'ModelPackageInputs.f90',
    'ModflowInput.f90',
    'SourceCommon.f90',
    'SourceLoad.F90',
    'mf6blockfile' / 'AsciiInputLoadType.f90',
    'mf6blockfile' / 'IdmMf6File.f90',
    'mf6blockfile' / 'LoadMf6File.f90',
    'mf6blockfile' / 'StressGridInput.f90',
    'mf6blockfile' / 'StressListInput.f90',
    'mf6blockfile' / 'StructArray.f90',
    'mf6blockfile' / 'StructVector.f90',
    # todo: include DFN selector files with generator instead?
    'selector' / 'IdmDfnSelector.f90',
    'selector' / 'IdmExgDfnSelector.f90',
    'selector' / 'IdmSimDfnSelector.f90',
)


foreach m : models
    mcc = m[0].to_upper() + m.substring(1, 3)  # camel case
    modflow_sources += files(
        'selector' / f'Idm@mcc@DfnSelector.f90',
    )
endforeach