[Mesh]
  type = FileMesh
  file = half_model_big.e
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 -9.81'
  displacements = 'disp_x disp_y disp_z'
  biot_coefficient = 1.0
  stress_free_temperature = 'temp_ref'
[]

[Variables]
  [temp]
    scaling = 1E-8
  []
  [pp]
    # scaling = 1E-5
  []
  [disp_x]
    scaling = 1E-12
  []
  [disp_y]
    scaling = 1E-12
  []
  [disp_z]
    scaling = 1E-12
  []
[]

[Functions]
  [dts]
    type = PiecewiseLinear
    x = '-2592000 -1 0 3600 86400 864000 864001 950400 1728000 1728001 1814400 2592000 2592001 2678400 3456000'
    y = '86400 86400 1 2400 43200   86400    1    43200 86400       1      43200   86400 1 43200 86400'
  []
  [pres_func]
    type = ParsedFunction
    expression = '(-(z-1000)) * rho * g' # z in mesh is negative and model is 1000 m below sea level
    symbol_names =  'rho g'
    symbol_values = '1000 9.81'
  []
  [temp_func]
      type = ParsedFunction
      expression = 't_surf + (-z) * 0.025'
      symbol_names =  't_surf'
      symbol_values = '323.15'
  []
  [lithostat_xx]
    type = ParsedFunction
    expression = '-(-(z-1000)) * 2500 * 9.81 * 0.9'
  []
  [lithostat_yy]
    type = ParsedFunction
    expression = '-(-(z-1000)) * 2500 * 9.81 * 1.1'
  []
  [lithostat_zz]
    type = ParsedFunction
    expression = '-(-(z-1000)) * 2500 * 9.81'
  []
  [lithostat_xx_force]
    type = ParsedFunction
    expression = '(-(z-1000)) * 2500 * 9.81 * 0.9'
  []
  [lithostat_yy_force]
    type = ParsedFunction
    expression = '(-(z-1000)) * 2500 * 9.81 * 1.1'
  []
  [lithostat_zz_force]
    type = ParsedFunction
    expression = '(-(z-1000)) * 2500 * 9.81'
  []
  [inject_hot]
    type = ParsedFunction
    expression = 'if(t >= 0 & t < 864000, 1, if(t >= 1728000 & t < 2592000, 1, 0))'
  []
  [produce_cold]
    type = ParsedFunction
    expression = 'if(t >= 0 & t < 864000, 1, if(t >= 1728000 & t < 2592000, 1, 0))'
  []
  [produce_hot]
    type = ParsedFunction
    expression = 'if(t >= 864000 & t < 1728000, 1, if(t >= 2592000 & t < 3456000, 1, 0))'
  []
  [inject_cold]
    type = ParsedFunction
    expression = 'if(t >= 864000 & t < 1728000, 1, if(t >= 2592000 & t < 3456000, 1, 0))'
  []
  [injection_rate_value]
    type = ParsedFunction
    symbol_names = true_screen_area
    symbol_values = true_screen_area
    expression = '-${inject_fluid_mass}/(true_screen_area * ${inject_time})'
  []
  [production_rate_value]
    type = ParsedFunction
    symbol_names = true_screen_area
    symbol_values = true_screen_area
    expression = '${produce_fluid_mass}/(true_screen_area * ${produce_time})'
  []
[]

inject_fluid_mass = 3456000
inject_time = 1728000
produce_fluid_mass = 3456000
produce_time = 1728000

[ICs]
  [pressure_ic]
    type = FunctionIC
    variable = pp
    function = pres_func
  []
  [temperature_ic]
    type = FunctionIC
    variable = temp
    function = temp_func
  []
[]

[BCs]
  [disp_x_null]
    type = DirichletBC
    variable = disp_x
    value = 0
    boundary = 'right_side'
  []
  [disp_y_null]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = 'front_side'
  []
  [disp_z_null]
    type = DirichletBC
    variable = disp_z
    value = 0
    boundary = 'bottom_side'
  []
  [stress_zz]
    type = Pressure
    boundary = 'top_side'
    variable = disp_z
    function = lithostat_zz_force
  []
  [stress_xx_max]
    type = Pressure
    boundary = 'left_side'
    variable = disp_x
    function = lithostat_xx_force
  []
  [stress_yy_max]
    type = Pressure
    boundary = 'back_side'
    variable = disp_y
    function = lithostat_yy_force
  []
  [pressure]
    type = FunctionDirichletBC
    variable = pp
    function = pres_func
    boundary = 'left_side back_side right_side'
  []
  [temperature]
    type = FunctionDirichletBC
    variable = temp
    function = temp_func
    boundary = 'top_side bottom_side'
  []
  [inject_heat]
    type = DirichletBC
    variable = temp
    boundary = 'hot_area'
    value = 373.15
  []
  [inject_fluid_hot]
    type = PorousFlowSink
    variable = pp
    boundary = 'hot_area'
    flux_function = injection_rate_value
  []
  [produce_heat]
    type = PorousFlowSink
    variable = temp
    boundary = 'hot_area'
    flux_function = production_rate_value
    fluid_phase = 0
    use_enthalpy = true
  []
  [produce_fluid_hot]
    type = PorousFlowSink
    variable = pp
    boundary = 'hot_area'
    flux_function = production_rate_value
  []
  [inject_cold]
    type = DirichletBC
    variable = temp
    boundary = 'cold_area'
    value = 303.15 # 30 °C
  []
  [inject_fluid_cold]
    type = PorousFlowSink
    variable = pp
    boundary = 'cold_area'
    flux_function = injection_rate_value
  []
  [produce_cold]
    type = PorousFlowSink
    variable = temp
    boundary = 'cold_area'
    flux_function = production_rate_value
    fluid_phase = 0
    use_enthalpy = true
  []
  [produce_fluid_cold]
    type = PorousFlowSink
    variable = pp
    boundary = 'cold_area'
    flux_function = production_rate_value
  []
[]

[Controls]
  [hot_inject_on]
    type = ConditionalFunctionEnableControl
    enable_objects = 'BCs::inject_heat BCs::inject_fluid_hot'
    conditional_function = inject_hot
    implicit = false
    execute_on = 'initial timestep_begin'
  []
  [hot_produce_on]
    type = ConditionalFunctionEnableControl
    enable_objects = 'BCs::produce_heat BCs::produce_fluid_hot'
    conditional_function = produce_hot
    implicit = false
    execute_on = 'initial timestep_begin'
  []
  [cold_inject_on]
    type = ConditionalFunctionEnableControl
    enable_objects = 'BCs::inject_cold BCs::inject_fluid_cold'
    conditional_function = inject_cold
    implicit = false
    execute_on = 'initial timestep_begin'
  []
  [cold_produce_on]
    type = ConditionalFunctionEnableControl
    enable_objects = 'BCs::produce_cold BCs::produce_fluid_cold'
    conditional_function = produce_cold
    implicit = false
    execute_on = 'initial timestep_begin'
  []
[]

[Kernels]
  [mass_dot]
    # type = PorousFlowMassTimeDerivative
    type = PorousFlowFullySaturatedMassTimeDerivative
    # fluid_component = 0
    variable = pp
  []
  [advection]
    type = PorousFlowFullySaturatedAdvectiveFlux
    # type = PorousFlowAdvectiveFlux # looks normal!
    fluid_component = 0
    variable = pp
    # biot_coefficient = 0.6
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = temp
  []
  [convection]
    type = PorousFlowFullySaturatedHeatAdvection
    variable = temp
  []
  [heat_conduction]
    type = PorousFlowHeatConduction
    variable = temp
  []
  [vol_strain_rate_water]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 0
    variable = pp
  []
  [vol_strain_rate_heat]
    type = PorousFlowHeatVolumetricExpansion
    variable = temp
  []
  [grad_stress_x]
    type = StressDivergenceTensors
    temperature = temp
    variable = disp_x
    eigenstrain_names = thermal_contribution
    use_displaced_mesh = false
    component = 0
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    use_displaced_mesh = false
    component = 0
  []
  [gravity]
    type = Gravity
    use_displaced_mesh = false
    variable = disp_z
    value = -10E-6 # MPa
  []
  [grad_stress_y]
    type = StressDivergenceTensors
    temperature = temp
    variable = disp_y
    eigenstrain_names = thermal_contribution
    use_displaced_mesh = false
    component = 1
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    use_displaced_mesh = false
    component = 1
  []
  [grad_stress_z]
    type = StressDivergenceTensors
    temperature = temp
    variable = disp_z
    eigenstrain_names = thermal_contribution
    use_displaced_mesh = false
    component = 2
  []
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    use_displaced_mesh = false
    component = 2
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'temp pp disp_x disp_y disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [produced_mass_water]
    type = PorousFlowSumQuantity
  []
  [produced_heat]
    type = PorousFlowSumQuantity
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 3E9
    viscosity = 1.0e-3
    density0 = 1000.0
    thermal_expansion = 7.5e-6
  []
[]

[Materials]
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 1E-10
    fluid_bulk_modulus = 3E9
  []
  [temperature]
    type = PorousFlowTemperature
    temperature = temp
  []
  [PS]
    type = PorousFlow1PhaseFullySaturated
    porepressure = pp
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [fp_mat]
    type = FluidPropertiesMaterialPT
    pressure = pp
    temperature = temp
    fp = simple_fluid
  []
  [porosity_cap]
    type = PorousFlowPorosityConst
    porosity = 0.05
    block = 'cap basement'
  []
  [porosity_res]
    type = PorousFlowPorosityConst
    porosity = 0.25
    block = 'res'
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 6.0E9
    poissons_ratio = 0.25
    # youngs_modulus = 1.0E10
  []
  [strain]
    type = ComputeSmallStrain
    eigenstrain_names = 'thermal_contribution initial_stress'
  []
  [initial_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'lithostat_xx 0 0  0 lithostat_yy 0  0 0 lithostat_zz'
    eigenstrain_name = initial_stress
  []
  [stress]
    type = ComputeLinearElasticStress # ComputeStrainIncrementBasedStress
  []
  [thermal_contribution]
    type = ComputeThermalExpansionEigenstrain
    temperature = temp
    thermal_expansion_coeff = 5E-6 # this is the linear thermal expansion coefficient
    eigenstrain_name = thermal_contribution
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
  []
  [rock_internal_energy]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 1200.0
    block = 'cap res basement'
  []
  [permeability_cap]
    type = PorousFlowPermeabilityConst
    permeability = '1E-17 0 0   0 1E-17 0   0 0 1E-17'
    block = 'cap basement'
  []
  [permeability_res]
    type = PorousFlowPermeabilityConst
    permeability = '4.9E-14 0 0   0 4.9E-14 0   0 0 4.9E-14' # 50 miliDarcy
    block = 'res'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
  []
  [thermal_conductivity_caps]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '1.6 0 0  0 1.6 0  0 0 1.6'
    block = 'cap basement'
  []
  [thermal_conductivity_res]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '2.2 0 0  0 2.2 0  0 0 2.2'
    block = 'res'
  []
[]

[Preconditioning]
  active = 'p1'
  [p1]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  []
  [mumps]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type'
    petsc_options_value = 'gmres      lu       mumps                         NONZERO'
  []
  [smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2'
  []
  [basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  []
  [preferred]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

[Executioner]
  type = Transient
  solve_type = Newton
  [TimeStepper]
    type = FunctionDT
    function = dts
  []
  # automatic_scaling = true
  # compute_scaling_once = false
  end_time = 3456000
  start_time = -2592000
  # nl_rel_step_tol = 1e-12
[]

[Outputs]
  exodus = true
  [csv]
    type = CSV
    interval = 1
  []
[]

[Postprocessors]
  [true_screen_area] # the injection/production area
    type = AreaPostprocessor
    boundary = hot_area
    execute_on = 'initial'
    outputs = 'none'
  []
  [hot_well_disp_x]
    type =  PointValue
    variable = disp_x
    point = '-999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [hot_well_disp_y]
    type =  PointValue
    variable = disp_y
    point = '-999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [hot_well_disp_z]
    type =  PointValue
    variable = disp_z
    point = '-999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [cold_well_disp_x]
    type =  PointValue
    variable = disp_x
    point = '999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [cold_well_disp_y]
    type =  PointValue
    variable = disp_y
    point = '999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [cold_well_disp_z]
    type =  PointValue
    variable = disp_z
    point = '999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [hot_well_press]
    type =  PointValue
    variable = pp
    point = '-999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [hot_well_temp]
    type = PointValue
    variable = temp
    point = '-999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [cold_well_press]
    type =  PointValue
    variable = pp
    point = '999 -1000 -51'
    execute_on = TIMESTEP_END
  []
  [cold_well_temp]
    type = PointValue
    variable = temp
    point = '999 -1000 -51'
    execute_on = TIMESTEP_END
  []
[]

[Debug]
  show_var_residual_norms = true
[]

[AuxVariables]
  [temp_ref]
  []
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [./p0]
    family = LAGRANGE
    order = FIRST
  [../]
  [./dP]
    family = LAGRANGE
    order = FIRST
  [../]
[]

[AuxKernels]
  [./p0_kernel]
    type = FunctionAux
    variable = 'p0'
    function = '-(1000*9.81*(z - 1000))'
    execute_on = 'initial'
  [../]
  [./dP_kernel]
    type = ParsedAux
    variable = 'dP'
    args = 'pp p0'
    function = 'pp-p0'
    execute_on = 'TIMESTEP_END'
  [../]
  [temp_ref]
    type = FunctionAux
    function = 'temp_func'
    variable = temp_ref
    execute_on = 'INITIAL'
  []
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  []
  [stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
  []
  [stress_yx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yx
    index_i = 1
    index_j = 0
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
  []
  [stress_zx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zx
    index_i = 2
    index_j = 0
  []
  [stress_zy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zy
    index_i = 2
    index_j = 1
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  []
[]
