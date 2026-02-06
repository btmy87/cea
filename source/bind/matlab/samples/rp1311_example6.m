%% rp1311_example6
% recreation of source\bind\c\samples\rp1311_example6.c
%------------------------------------------------------------------
% Problem Specification
%------------------------------------------------------------------

% Species
reactants       = ["H2(L)", "O2(L)"];
fuel_weights    = [    1.0,     0.0];
oxidant_weights = [    0.0,     1.0];
nr = length(reactants);

% Products
omitted_products = strings();

% Thermo States
p0 = 1.0;
T0 = 298.15;
eq_ratio = 1.0;


%----------------------------------------------------------------
% CEA Setup
%----------------------------------------------------------------

cea = CEA();
cea.set_log_level("CEA_LOG_CRITICAL");

% Mixtures
reac = cea.mixture_create(reactants);
prod = cea.mixture_create_from_reactants( ...
    reactants, omitted_products );

% Solver
solver = cea.detonation_solver_create_with_reactants(prod, reac);

% Solution
soln = cea.detonation_solution_create();


%----------------------------------------------------------------
% Solve the detonation problem
%----------------------------------------------------------------

of_ratio = cea.mixture_chem_eq_ratio_to_of_ratio(reac, nr, oxidant_weights, fuel_weights, eq_ratio);
weights = cea.mixture_of_ratio_to_weights(reac, nr, oxidant_weights, fuel_weights, of_ratio);

cea.detonation_solver_solve(solver, soln, weights, T0, p0, false);

temperature    = cea.detonation_solution_get_property(soln, "CEA_DETONATION_TEMPERATURE"   , 1);
pressure       = cea.detonation_solution_get_property(soln, "CEA_DETONATION_PRESSURE"      , 1);
velocity       = cea.detonation_solution_get_property(soln, "CEA_DETONATION_VELOCITY"      , 1);
mach           = cea.detonation_solution_get_property(soln, "CEA_DETONATION_MACH"          , 1);
sonic_velocity = cea.detonation_solution_get_property(soln, "CEA_DETONATION_SONIC_VELOCITY", 1);
gamma          = cea.detonation_solution_get_property(soln, "CEA_DETONATION_GAMMA"         , 1);
enthalpy       = cea.detonation_solution_get_property(soln, "CEA_DETONATION_ENTHALPY"      , 1);

fprintf( ...
    "%10s %10s %12s %8s %12s %10s %8s\n", ...
    "T (K)", "P (bar)", "Velocity (m/s)", "Mach", "Sonic Velocity (m/s)", "Gamma", "H/R" ...
);

fprintf( ...
    "%10.3f  %9.3f  %13.3f  %7.4f  %19.3f  %9.4f  %6.3f\n", ...
    temperature, pressure, velocity, mach, sonic_velocity, gamma, enthalpy ...
);

%----------------------------------------------------------------
% CEA Cleanup
%----------------------------------------------------------------
cea.detonation_solution_destroy(soln);
cea.detonation_solver_destroy(solver);
cea.mixture_destroy(prod);
cea.mixture_destroy(reac);
