#include "cea.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>


#define LEN(x) (sizeof(x) / sizeof((x)[0]))
#define ATM 1.01325

static int compare_close(cea_real a, cea_real b) {
  cea_real diff = a - b;
  cea_real denom = fmax(1.0, fabs(b));
  return fabs(diff) <= 1e-8 * denom;
}

#define CHECK_PROP(name, get_enum) \
  do { \
    cea_real val_direct = soln->name; \
    cea_real val_get; \
    cea_eqsolution_get_property(soln, get_enum, &val_get); \
    if (!compare_close(val_direct, val_get)) { \
      printf("%s MISMATCH: direct=%e getter=%e\n", #name, val_direct, val_get); \
      overall_fail = 1; \
    } \
  } while (0)

int main(void) {
  //------------------------------------------------------------------
  // TP Problem Specification
  //------------------------------------------------------------------

  // Reactants
  const cea_string reactants[] = {"H2", "Air"};
  const cea_real fuel_moles[] = {1.0, 0.0};
  const cea_real oxidant_moles[] = {0.0, 1.0};

  // Products
  const cea_string products[] = {
      "Ar",   "C", "CO", "CO2", "H",  "H2",   "H2O", "HNO", "HO2", "HNO2",
      "HNO3", "N", "NH", "NO",  "N2", "N2O3", "O",   "O2",  "OH",  "O3"};

  // Mixture States
  const cea_real pressures[] = {1.00 * ATM};
  const cea_real temperatures[] = {3000.0};
  const cea_real chem_eq_ratios[] = {1.0};

  //------------------------------------------------------------------
  // CEA Setup
  //------------------------------------------------------------------

  cea_set_log_level(CEA_LOG_ERROR);
  cea_init();

  // Mixtures
  cea_mixture reac, prod;
  cea_mixture_create(&reac, LEN(reactants), reactants);
  cea_mixture_create(&prod, LEN(products), products);

  // EqSolver
  cea_eqsolver solver;
  cea_solver_opts opts;
  cea_solver_opts_init(&opts);
  opts.reactants = reac;
  cea_eqsolver_create_with_options(&solver, prod, opts);

  // EqSolution
  cea_eqsolution soln;
  cea_eqsolution_create(&soln, solver);

  // EqPartials
  cea_eqpartials partials;
  cea_eqpartials_create(&partials, solver);

  //------------------------------------------------------------------
  // Unit Conversions
  //------------------------------------------------------------------

  cea_real fuel_weights[LEN(reactants)];
  cea_mixture_moles_to_weights(reac, LEN(reactants), fuel_moles, fuel_weights);
  cea_real oxidant_weights[LEN(reactants)];
  cea_mixture_moles_to_weights(reac, LEN(reactants), oxidant_moles,
                               oxidant_weights);

  cea_real of_ratios[LEN(chem_eq_ratios)];
  for (int ir = 0; ir < LEN(chem_eq_ratios); ++ir) {
    cea_mixture_chem_eq_ratio_to_of_ratio(reac, LEN(reactants), oxidant_weights,
                                          fuel_weights, chem_eq_ratios[ir],
                                          &of_ratios[ir]);
  }

  //------------------------------------------------------------------
  // Equilibrium Solve
  //------------------------------------------------------------------

  int overall_fail = 0;

  for (int ir = 0; ir < LEN(of_ratios); ++ir) {
    cea_real weights[LEN(reactants)];
    cea_mixture_of_ratio_to_weights(reac, LEN(reactants), oxidant_weights,
                                    fuel_weights, of_ratios[ir], weights);

    for (int ip = 0; ip < LEN(pressures); ++ip) {
      for (int it = 0; it < LEN(temperatures); ++it) {
        cea_eqsolver_solve_with_partials(solver, CEA_TP, temperatures[it],
                                         pressures[ip], weights, soln,
                                         partials);

        CHECK_PROP(T, CEA_TEMPERATURE);
        CHECK_PROP(pressure, CEA_PRESSURE);
        CHECK_PROP(volume, CEA_VOLUME);
        CHECK_PROP(density, CEA_DENSITY);
        CHECK_PROP(M, CEA_M);
        CHECK_PROP(MW, CEA_MW);
        CHECK_PROP(enthalpy, CEA_ENTHALPY);
        CHECK_PROP(energy, CEA_ENERGY);
        CHECK_PROP(entropy, CEA_ENTROPY);
        CHECK_PROP(gibbs_energy, CEA_GIBBS_ENERGY);
        CHECK_PROP(gamma_s, CEA_GAMMA_S);
        CHECK_PROP(cp_fr, CEA_FROZEN_CP);
        CHECK_PROP(cv_fr, CEA_FROZEN_CV);
        CHECK_PROP(cp_eq, CEA_EQUILIBRIUM_CP);
        CHECK_PROP(cv_eq, CEA_EQUILIBRIUM_CV);
        CHECK_PROP(viscosity, CEA_VISCOSITY);
        CHECK_PROP(conductivity_fr, CEA_FROZEN_CONDUCTIVITY);
        CHECK_PROP(conductivity_eq, CEA_EQUILIBRIUM_CONDUCTIVITY);
        CHECK_PROP(Pr_fr, CEA_FROZEN_PRANDTL);
        CHECK_PROP(Pr_eq, CEA_EQUILIBRIUM_PRANDTL);
      }
    }
  }

  //----------------------------------------------------------------
  // CEA Cleanup
  //----------------------------------------------------------------
  cea_eqpartials_destroy(&partials);
  cea_eqsolution_destroy(&soln);
  cea_eqsolver_destroy(&solver);
  cea_mixture_destroy(&prod);
  cea_mixture_destroy(&reac);

  return overall_fail;
}
