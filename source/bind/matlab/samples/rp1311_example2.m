%% rp1311_example2.m
% recreate example from source/bind/c/rp1311_example2.m
close all
clear
clc

ATM = 1.01325;

%% Define Problem
reactants = ["H2", "Air"];
fuel_weights = [1.0, 0.0];
oxidant_weights = [0.0, 1.0];
products = ["Ar",    "C",    "CO",   "CO2",  "H", ...
            "H2",    "H2O",  "HNO",  "HO2",  "HNO2", ...
            "HNO3",  "N",    "NH",   "NO",   "N2", ...
            "N2O3",  "O",    "O2",   "OH",   "O3"];

densities = [1.0e3*9.1864e-5, 1.0e3*8.0877e-6, 1.0e3*6.6054e-7 ];
temperatures = 3000;
phi_eq_ratios = 1.0;

%% cea setup
cea = CEA();
cea.set_log_level("CEA_LOG_CRITICAL"); % lack of a console causes a crash if it actually writes a log
                                       % LOG_NONE doesn't do what you'd expect

% mixtures
reac = cea.mixture_create(reactants);
prod = cea.mixture_create(products);

% eqsolver
opts = cea.solver_opts_init();
opts.Value.reactants = reac;
opts.Value.transport = 1; % this isn't taking?

solver = cea.eqsolver_create_with_options(prod, opts);

% eqsolution
soln = cea.eqsolution_create(solver);

% partials
partials = cea.eqpartials_create(solver);

%% unit conversions
of_ratios = zeros(size(phi_eq_ratios));
for i = 1:length(of_ratios)
    of_ratios(i) = cea.mixture_chem_eq_ratio_to_of_ratio(reac, ...
        length(reactants), oxidant_weights, fuel_weights, ...
        phi_eq_ratios(i));
end

%% equilibrium solve
fprintf("%10s  %10s  %10s  %12s  %12s  %12s\n",...
        "O/F Ratio", "P (atm)", "T (K)",...
        "H (cal/g)", "Cp (cal/g-K)", "Viscosity (mP)");

for ir = 1:length(of_ratios)
    weights = cea.mixture_of_ratio_to_weights(reac, length(reactants), ...
        oxidant_weights, fuel_weights, of_ratios(ir));
    for ip = 1:length(densities)
        for it = 1:length(temperatures)
            cea.eqsolver_solve_with_partials(solver, "CEA_TV", ...
                temperatures(it), 1.0./densities(ip), ...
                weights, soln, partials);
            pressure = cea.eqsolution_get_property(soln, "CEA_PRESSURE");
            enthalpy = cea.eqsolution_get_property(soln, "CEA_ENTHALPY");
            heat_capacity = cea.eqsolution_get_property(soln, "CEA_EQUILIBRIUM_CP");
            viscosity = cea.eqsolution_get_property(soln, "CEA_VISCOSITY");
            pressure = pressure./ ATM;
            enthalpy = enthalpy./4.184;
            heat_capacity = heat_capacity./4.184;
            fprintf("%10.2f  %10.2f  %10.2f  %12.5e  %10.5f  %10.5f\n", ...
                of_ratios(ir), pressure, temperatures(it), ...
                enthalpy, heat_capacity, viscosity);

        end
    end
end

%% cleanup
cea.eqpartials_destroy(partials);
cea.eqsolution_destroy(soln);
cea.eqsolver_destroy(solver);
cea.mixture_destroy(prod);
cea.mixture_destroy(reac);