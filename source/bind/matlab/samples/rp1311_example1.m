%% rp1311_example1.m
% recreate rp1311_example1.c
close all
clear
clc

%% Load CEA
load_cea();
cea = CEA();
cea.set_log_level("CEA_LOG_CRITICAL"); % lack of a console causes a crash if it actually writes a log
                                       % LOG_NONE doesn't do what you'd expect

%% Define Problem
reactants = ["H2", "Air"];
fuel_moles = [1.0, 0.0];
oxidant_moles = [0.0, 1.0];
products = ["Ar",    "C",    "CO",   "CO2",  "H", ...
            "H2",    "H2O",  "HNO",  "HO2",  "HNO2", ...
            "HNO3",  "N",    "NH",   "NO",   "N2", ...
            "N2O3",  "O",    "O2",   "OH",   "O3"];

pressures = [1.0, 0.1, 0.01].*1.01325;
temperatures = [3000.0, 2000.0];
chem_eq_ratios = [1.0, 1.5];

% create mixtures
reac = cea.mixture_create(reactants);
prod = cea.mixture_create(products);

% create solver
opts = cea.solver_opts_init();
opts.Value.reactants = reac;
solver = cea.eqsolver_create_with_options(prod, opts);

% equilibrium solution
soln = cea.eqsolution_create(solver);

% partials
partials = cea.eqpartials_create(solver);

% unit conversions
fuel_weights = cea.mixture_moles_to_weights(reac, ...
    length(reactants), fuel_moles);
oxidant_weights = cea.mixture_moles_to_weights(reac, ...
    length(reactants), oxidant_moles);

of_ratios = zeros(size(chem_eq_ratios));
for i = 1:length(of_ratios)
    of_ratios(i) = cea.mixture_chem_eq_ratio_to_of_ratio( ...
        reac, length(reactants), oxidant_weights, fuel_weights, ...
        chem_eq_ratios(i));
end

% equilibrium solve
fprintf("%10s  %10s  %10s  %10s  %12s  %12s  %12s  %12s\n", ...
        "T (K)", "P (Pa)", "Chem Equiv", "O/F Ratio", ...
        "H2 (wtfrac)", "Air (wtfrac)", "H (cal/g)", "Cp (cal/g-K)");

for ir = 1:length(of_ratios)
    weights = cea.mixture_of_ratio_to_weights(reac, length(reactants), ...
        oxidant_weights, fuel_weights, of_ratios(ir));

    for ip = 1:length(pressures)
        for it = 1:length(temperatures)
            cea.eqsolver_solve_with_partials(solver, "CEA_TP", ...
                temperatures(it), pressures(ip), weights, soln, partials);
            h = cea.eqsolution_get_property(soln, "CEA_ENTHALPY");
            cp = cea.eqsolution_get_property(soln, "CEA_EQUILIBRIUM_CP");
            h = h/4.184;
            cp = cp/4.184;

            fprintf("%10.2f  %10.2f  %10.2f  %10.6f  %12.5e  %12.5e  %12.5e  %12.5e\n", ...
                temperatures(it), pressures(ip), chem_eq_ratios(ir), ...
                of_ratios(ir), weights(1), weights(2), h, cp);
        end
    end
end

cea.eqpartials_destroy(partials);
%cea.eqsolution_destroy(soln);
cea.eqsolver_destroy(solver);
cea.mixture_destroy(prod);
cea.mixture_destroy(reac);


