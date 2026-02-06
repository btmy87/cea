BAR = 10000.00;
R = 8314.51;

%------------------------------------------------------------------
% Problem Specification
%------------------------------------------------------------------

% Species
reactants       = ["H2(L)", "O2(L)"];
reactant_temps  = [  20.27,   90.17];
fuel_weights    = [    1.0,     0.0];
oxidant_weights = [    0.0,     1.0];
nr = length(reactants);

% Products
omitted_products = strings([]);

% Thermo States
pressures = 53.3172*BAR;
of_ratio  = 5.55157;

%----------------------------------------------------------------
% CEA Setup
%----------------------------------------------------------------

cea = CEA();
cea.set_log_level("CEA_LOG_CRITICAL");

% Mixtures
reac = cea.mixture_create(reactants);
prod = cea.mixture_create_from_reactants(reactants, omitted_products);

% Solver
solver = cea.rocket_solver_create_with_reactants(prod, reac);

% Solution
num_pts = 5;
soln = cea.rocket_solution_create(solver);

%----------------------------------------------------------------
% Rocket Solve
%----------------------------------------------------------------
weights = cea.mixture_of_ratio_to_weights(reac, nr, oxidant_weights, fuel_weights, of_ratio);

pip   = [ 10.0 ];
subar = [ 1.58 ];
supar = [ 25.0 ];
pc = 53.3172;

hc = cea.mixture_calc_property_multitemp(reac, "CEA_ENTHALPY", weights, reactant_temps);
hc = hc/R;

cea.rocket_solver_solve_iac(solver, soln, weights, pc, pip, 1, subar, 1, supar, 1, 0, hc, true, 0.0, false);

temperature = cea.rocket_solution_get_property(soln, "CEA_ROCKET_TEMPERATURE"   , num_pts);
pressure    = cea.rocket_solution_get_property(soln, "CEA_ROCKET_PRESSURE"      , num_pts);
gamma       = cea.rocket_solution_get_property(soln, "CEA_ROCKET_GAMMA_S"       , num_pts);
mw          = cea.rocket_solution_get_property(soln, "CEA_ROCKET_MW"            , num_pts);
mach        = cea.rocket_solution_get_property(soln, "CEA_MACH"                 , num_pts);
area_ratio  = cea.rocket_solution_get_property(soln, "CEA_AE_AT"                , num_pts);
isp         = cea.rocket_solution_get_property(soln, "CEA_ISP"                  , num_pts);
cstar       = cea.rocket_solution_get_property(soln, "CEA_C_STAR"               , num_pts);
cf          = cea.rocket_solution_get_property(soln, "CEA_COEFFICIENT_OF_THRUST", num_pts);

fprintf( ...
    "%10s %10s %10s %10s %10s %10s %10s %10s %10s \n", ...
    "T (K)", "P (bar)", "gamma_s", "MW", "Mach", "Ae/At", "ISP", "C* (m/s)", "C_f" ...
);

for i = 1:num_pts
    fprintf( ...
        "%10.4f  %10.4f  %10.4f  %10.4f %10.4f  %10.4e  %10.4e  %10.4e  %10.4e\n", ...
        temperature(i), pressure(i), gamma(i), mw(i), mach(i), ...
        area_ratio(i), isp(i), cstar(i), cf(i) ...
    );
end

%----------------------------------------------------------------
% CEA Cleanup
%----------------------------------------------------------------
cea.rocket_solution_destroy(soln);
cea.rocket_solver_destroy(solver);
cea.mixture_destroy(prod);
cea.mixture_destroy(reac);
