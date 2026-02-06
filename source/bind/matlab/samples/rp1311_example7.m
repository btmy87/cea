%------------------------------------------------------------------
% Problem Specification
%------------------------------------------------------------------

% Species
reactants     = ["H2",  "O2",  "Ar" ];
moles         = [0.050, 0.050, 0.900];
nr = length(reactants);

% Products
omitted_products = strings();

% Thermo States
p0 = 0.1;
T0 = 300.0;

% Shock Parameters
u1 = 1400.0;


%----------------------------------------------------------------
% CEA Setup
%----------------------------------------------------------------

cea = CEA();
cea.set_log_level("CEA_LOG_CRITICAL");

% Mixtures
reac = cea.mixture_create(reactants);
prod = cea.mixture_create_from_reactants(reactants, omitted_products);

% Solver
solver = cea.shock_solver_create_with_reactants(prod, reac);

% Solution
num_pts = 3;
soln = cea.shock_solution_create(num_pts);


%----------------------------------------------------------------
% Solve the shock problem
%----------------------------------------------------------------
weights = cea.mixture_moles_to_weights(reac, nr, moles);

cea.shock_solver_solve(solver, soln, weights, T0, p0, u1, false, true, false, false);

temperature = cea.shock_solution_get_property(soln, "CEA_SHOCK_TEMPERATURE"   , num_pts);
pressure    = cea.shock_solution_get_property(soln, "CEA_SHOCK_PRESSURE"      , num_pts);
velocity    = cea.shock_solution_get_property(soln, "CEA_SHOCK_VELOCITY"      , num_pts);
mach        = cea.shock_solution_get_property(soln, "CEA_SHOCK_MACH"          , num_pts);
v_sonic     = cea.shock_solution_get_property(soln, "CEA_SHOCK_SONIC_VELOCITY", num_pts);

fprintf( ...
    "%14s %14s %14s %14s %14s \n", ...
    "T (K)", "P (bar)", "u (m/s)", "Mach", "Son. Vel. (m/s)" ...
);

for i = 1:num_pts
    fprintf( ...
        "%14.4f  %14.4f  %14.4f  %14.4f  %14.4f\n", ...
        temperature(i), pressure(i), velocity(i), mach(i), v_sonic(i) ...
    );
end

%----------------------------------------------------------------
% CEA Cleanup
%----------------------------------------------------------------
cea.shock_solution_destroy(soln);
cea.shock_solver_destroy(solver);
cea.mixture_destroy(prod);
cea.mixture_destroy(reac);

