#include "stdbool.h"
#include "cea_enum.h"

#ifdef __cplusplus
extern "C" {
#endif

// Enumerations
typedef enum { CEA_LOG_LEVEL_ENUM                } cea_log_level;
typedef enum { CEA_EQUILIBRIUM_TYPE_ENUM         } cea_equilibrium_type;
typedef enum { CEA_EQUILIBRIUM_SIZE_ENUM         } cea_equilibrium_size;
typedef enum { CEA_PROPERTY_TYPE_ENUM            } cea_property_type;
typedef enum { CEA_ROCKET_PROPERTY_TYPE_ENUM     } cea_rocket_property_type;
typedef enum { CEA_SHOCK_PROPERTY_TYPE_ENUM      } cea_shock_property_type;
typedef enum { CEA_DETONATION_PROPERTY_TYPE_ENUM } cea_detonation_property_type;
typedef enum { CEA_ERROR_CODE_ENUM               } cea_error_code;

// Type Shorthand, Opaque Types as void*
typedef void*              cea_reactant;
typedef void*              cea_mixture;
typedef void*              cea_eqsolver;
//typedef void*              cea_eqsolution;
typedef void*              cea_eqpartials;
typedef void*              cea_rocket_solver;
typedef void*              cea_rocket_solution;
typedef void*              cea_shock_solver;
typedef void*              cea_shock_solution;
typedef void*              cea_detonation_solver;
typedef void*              cea_detonation_solution;
typedef cea_error_code     cea_err;
typedef const char*        cea_string;
typedef int                cea_int;
typedef double             cea_real;
typedef double*            cea_array;

typedef struct {
    cea_string name;
    cea_int num_elements;
    const cea_string *elements;
    const cea_real *coefficients;
    char has_molecular_weight;
    cea_real molecular_weight;
    char has_enthalpy;
    cea_real enthalpy;
    cea_string enthalpy_units;
    char has_temperature;
    cea_real temperature;
    cea_string temperature_units;
} cea_reactant_input;

// fortran array descriptor format
// int64 addr
// int64 elem_size
// int64 reserved
// int64 flags
// int64 dimensions
// int64 reserved
// int64 num_elems, dim 1
// int64 mem_dist, dim 1
// int64, lower_bound, dim 1
// last 3 fields are repeated for each dimension
typedef long long int64;
typedef struct {
        double T;
        double *nj;
        int64 nj_elem_size;
        int64 nj_reserved1;
        int64 nj_flags;
        int64 nj_dimensions;
        int64 nj_reserved2;
        int64 nj_num1;
        int64 nj_dist1;
        int64 nj_lb1;
        double n;
        int thermo_num_species; int thermo_num_species_pad;
        double *thermo_cp;
        int64 thermo_cp_elem_size;
        int64 thermo_cp_reserved1;
        int64 thermo_cp_flags;
        int64 thermo_cp_dimensions;
        int64 thermo_cp_reserved2;
        int64 thermo_cp_num1;
        int64 thermo_cp_dist1;
        int64 thermo_cp_lb1;
        double *thermo_cv;
        int64 thermo_cv_elem_size;
        int64 thermo_cv_reserved1;
        int64 thermo_cv_flags;
        int64 thermo_cv_dimensions;
        int64 thermo_cv_reserved2;
        int64 thermo_cv_num1;
        int64 thermo_cv_dist1;
        int64 thermo_cv_lb1;
        double *thermo_enthalpy;
        int64 thermo_enthalpy_elem_size;
        int64 thermo_enthalpy_reserved1;
        int64 thermo_enthalpy_flags;
        int64 thermo_enthalpy_dimensions;
        int64 thermo_enthalpy_reserved2;
        int64 thermo_enthalpy_num1;
        int64 thermo_enthalpy_dist1;
        int64 thermo_enthalpy_lb1;
        double *thermo_entropy;
        int64 thermo_entropy_elem_size;
        int64 thermo_entropy_reserved1;
        int64 thermo_entropy_flags;
        int64 thermo_entropy_dimensions;
        int64 thermo_entropy_reserved2;
        int64 thermo_entropy_num1;
        int64 thermo_entropy_dist1;
        int64 thermo_entropy_lb1;
        double *thermo_energy;
        int64 thermo_energy_elem_size;
        int64 thermo_energy_reserved1;
        int64 thermo_energy_flags;
        int64 thermo_energy_dimensions;
        int64 thermo_energy_reserved2;
        int64 thermo_energy_num1;
        int64 thermo_energy_dist1;
        int64 thermo_energy_lb1;
        double *ln_nj;
        int64 ln_nj_elem_size;
        int64 ln_nj_reserved1;
        int64 ln_nj_flags;
        int64 ln_nj_dimensions;
        int64 ln_nj_reserved2;
        int64 ln_nj_num1;
        int64 ln_nj_dist1;
        int64 ln_nj_lb1;
        double *pi;
        int64 pi_elem_size;
        int64 pi_reserved1;
        int64 pi_flags;
        int64 pi_dimensions;
        int64 pi_reserved2;
        int64 pi_num1;
        int64 pi_dist1;
        int64 pi_lb1;
        double *pi_prev;
        int64 pi_prev_elem_size;
        int64 pi_prev_reserved1;
        int64 pi_prev_flags;
        int64 pi_prev_dimensions;
        int64 pi_prev_reserved2;
        int64 pi_prev_num1;
        int64 pi_prev_dist1;
        int64 pi_prev_lb1;
        double pi_e;
        double dpi_e;
        double *dln_nj;
        int64 dln_nj_elem_size;
        int64 dln_nj_reserved1;
        int64 dln_nj_flags;
        int64 dln_nj_dimensions;
        int64 dln_nj_reserved2;
        int64 dln_nj_num1;
        int64 dln_nj_dist1;
        int64 dln_nj_lb1;
        double *dnj_c;
        int64 dnj_c_elem_size;
        int64 dnj_c_reserved1;
        int64 dnj_c_flags;
        int64 dnj_c_dimensions;
        int64 dnj_c_reserved2;
        int64 dnj_c_num1;
        int64 dnj_c_dist1;
        int64 dnj_c_lb1;
        double dln_n;
        double dln_T;
        double *G;
        int64 G_elem_size;
        int64 G_reserved1;
        int64 G_flags;
        int64 G_dimensions;
        int64 G_reserved2;
        int64 G_num1;
        int64 G_dist1;
        int64 G_lb1;
        int64 G_num2;
        int64 G_dist2;
        int64 G_lb2;
        char type[8];
        double state1;
        double state2;
        double *b0;
        int64 b0_elem_size;
        int64 b0_reserved1;
        int64 b0_flags;
        int64 b0_dimensions;
        int64 b0_reserved2;
        int64 b0_num1;
        int64 b0_dist1;
        int64 b0_lb1;
        int *is_active;
        int64 is_active_elem_size;
        int64 is_active_reserved1;
        int64 is_active_flags;
        int64 is_active_dimensions;
        int64 is_active_reserved2;
        int64 is_active_num1;
        int64 is_active_dist1;
        int64 is_active_lb1;
        int j_liq; int j_sol;
        int j_switch; int last_cond_idx;
        int gas_converged; int condensed_converged;
        int moles_converged; int element_converged;
        int temperature_converged; int entropy_converged;
        int pi_converged; int ions_converged;
        int converged; int times_converged;
        double *mole_fractions;
        int64 mole_fractions_elem_size;
        int64 mole_fractions_reserved1;
        int64 mole_fractions_flags;
        int64 mole_fractions_dimensions;
        int64 mole_fractions_reserved2;
        int64 mole_fractions_num1;
        int64 mole_fractions_dist1;
        int64 mole_fractions_lb1;
        double *mass_fractions;
        int64 mass_fractions_elem_size;
        int64 mass_fractions_reserved1;
        int64 mass_fractions_flags;
        int64 mass_fractions_dimensions;
        int64 mass_fractions_reserved2;
        int64 mass_fractions_num1;
        int64 mass_fractions_dist1;
        int64 mass_fractions_lb1;
        double density; // (kg/m^3)
        double pressure; // bar
        double volume; // specific volume
        double M; // 1/n, Eq. 2.3a/b
        double MW; // Eq 2.4 a/b
        double enthalpy; // kJ/kg
        double energy;
        double gibbs_energy;
        double entropy; // kJ/kg-K
        double gamma_s;
        double viscosity;
        double cp_fr;
        double cp_eq;
        double cv_fr;
        double cv_eq;
        double conductivity_fr;
        double conductivity_eq;
        double Pr_fr;
        double Pr_eq;
} cea_eqsolution_t;
typedef cea_eqsolution_t* cea_eqsolution;

// Struct types for optional arguments
typedef struct {
    cea_real trace;
    char ions;
    char transport;
    cea_mixture reactants;
    cea_int ninsert;
    const cea_string* insert;
} cea_solver_opts;

// Initialize optional arguments
cea_err cea_solver_opts_init(cea_solver_opts *opts);

// Input string ownership
// All input strings/arrays (species, reactants, omit, insert, reactant_input fields) are copied by the
// C/Fortran layer during the call. Callers retain ownership and may free or reuse their buffers after
// the call returns.

// Species name utilities
cea_err cea_species_name_len(cea_int *name_len);

cea_int cea_test_add(cea_int a, cea_int b);

// Version Information
cea_err cea_version_major(cea_int *major);
cea_err cea_version_minor(cea_int *minor);
cea_err cea_version_patch(cea_int *patch);

// Logging Control
cea_err cea_set_log_level(const cea_log_level level);

// Initialization (not thread safe)
cea_err cea_init();
cea_err cea_init_thermo(const cea_string thermofile);
cea_err cea_init_trans(const cea_string transfile);
cea_err cea_is_initialized(cea_int *initialized);

//----------------------------------------------------------------------
// Mixture API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_mixture_create(
    cea_mixture *mix,
    const cea_int nspecies,
    const char** species
);

cea_err cea_mixture_create_w_ions(
    cea_mixture *mix,
    const cea_int nspecies,
    const char** species
);

cea_err cea_mixture_create_from_reactants(
    cea_mixture *mix,
    const cea_int nreactants,
    const char** reactants,
    const cea_int nomit,
    const char** omit
);

cea_err cea_mixture_create_from_reactants_w_ions(
    cea_mixture *mix,
    const cea_int nreactants,
    const char** reactants,
    const cea_int nomit,
    const char** omit
);

cea_err cea_mixture_create_from_input_reactants(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[]
);

cea_err cea_mixture_create_from_input_reactants_w_ions(
    cea_mixture *mix,
    const cea_int nreactants,
    const cea_reactant_input reactants[]
);

cea_err cea_mixture_destroy(
    cea_mixture *mix
);

// Get values
cea_err cea_mixture_get_num_species(
    const cea_mixture mix,
    cea_int *num_species
);

// Deprecated: returns a heap-allocated string that must be freed with cea_string_free.
cea_err cea_mixture_get_species_name(
    const cea_mixture *mix,
    const cea_int i_species,
    cea_string *species
);

// Deprecated: returns heap-allocated strings that must be freed with cea_string_array_free.
cea_err cea_mixture_get_species_names(
    const cea_mixture *mix,
    const cea_int nspecies,
    cea_string *species[]
);

// Buffer-based species name retrieval
cea_err cea_mixture_get_species_name_buf(
    const cea_mixture *mix,
    const cea_int i_species,
    char *species,
    const cea_int buf_len
);

cea_err cea_mixture_get_species_names_buf(
    const cea_mixture *mix,
    const cea_int nspecies,
    char *species,
    const cea_int stride
);

// Free functions for deprecated allocating getters
cea_err cea_string_free(cea_string species);
cea_err cea_string_array_free(cea_string species[], const cea_int nspecies);

// Unit Conversion
cea_err cea_mixture_moles_to_weights(
    const cea_mixture mix,
    const cea_int len,
    const cea_real moles[],
    cea_real weights[]
);

cea_err cea_mixture_weights_to_moles(
    const cea_mixture mix,
    const cea_int len,
    const cea_real weights[],
    cea_real moles[]
);

cea_err cea_mixture_per_mole_to_per_weight(
    const cea_mixture mix,
    const cea_int len,
    const cea_real per_mole[],
    cea_real per_weight[]
);

cea_err cea_mixture_per_weight_to_per_mole(
    const cea_mixture mix,
    const cea_int len,
    const cea_real per_weight[],
    cea_real per_mole[]
);

// Fuel/Oxidant Mixture Conversion
cea_err cea_mixture_chem_eq_ratio_to_of_ratio(
    const cea_mixture mix,
    const cea_int len,
    const cea_real oxidant_weights[],
    const cea_real fuel_weights[],
    const cea_real chem_eq_ratio,
    cea_real *of_ratio
);

cea_err cea_mixture_weight_eq_ratio_to_of_ratio(
    const cea_mixture mix,
    const cea_int len,
    const cea_real oxidant_weights[],
    const cea_real fuel_weights[],
    const cea_real weight_eq_ratio,
    cea_real *of_ratio
);

cea_err cea_mixture_of_ratio_to_weights(
    const cea_mixture mix,
    const cea_int len,
    const cea_real oxidant_weights[],
    const cea_real fuel_weights[],
    const cea_real of_ratio,
    cea_real reactant_weights[]
);

// Property Calculation
cea_err cea_mixture_calc_property(
    const cea_mixture mix,
    const cea_property_type type,
    const cea_int len_weights,
    const cea_real weights[],
    const cea_real temperature,
    cea_real *value
);

cea_err cea_mixture_calc_property_multitemp(
    const cea_mixture mix,
    const cea_property_type type,
    const cea_int len_weights,
    const cea_real weights[],
    const cea_int len_temperatures,
    const cea_real temperatures[],
    cea_real *value
);

cea_err cea_mixture_calc_property_tp(
    const cea_mixture mix,
    const cea_property_type type,
    const cea_int len_weights,
    const cea_real weights[],
    const cea_real temperature,
    const cea_real pressure,
    cea_real *value
);

cea_err cea_mixture_calc_property_tp_multitemp(
    const cea_mixture mix,
    const cea_property_type type,
    const cea_int len_weights,
    const cea_real weights[],
    const cea_int len_temperatures,
    const cea_real temperatures[],
    const cea_real pressure,
    cea_real *value
);


//----------------------------------------------------------------------
// Equilibrium Solver API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_eqsolver_create(
    cea_eqsolver *solver,
    const cea_mixture products
);

cea_err cea_eqsolver_create_with_reactants(
    cea_eqsolver *solver,
    const cea_mixture products,
    const cea_mixture reactants
);

cea_err cea_eqsolver_create_with_options(
    cea_eqsolver *solver,
    const cea_mixture products,
    const cea_solver_opts options
);

cea_err cea_eqsolver_destroy(
    cea_eqsolver *solver
);

// Solve
cea_err cea_eqsolver_solve(
    const cea_eqsolver solver,
    const cea_equilibrium_type type,
    const cea_real state1,
    const cea_real state2,
    const cea_array amounts,
    cea_eqsolution soln
);

cea_err cea_eqsolver_solve_with_partials(
    const cea_eqsolver solver,
    const cea_equilibrium_type type,
    const cea_real state1,
    const cea_real state2,
    const cea_array amounts,
    cea_eqsolution soln,
    cea_eqpartials eqpartials
);

// Querry functions
cea_err cea_eqsolver_get_size(
    const cea_eqsolver solver,
    const cea_equilibrium_size eq_variable,
    cea_int *value
);

// // Function interface
// cea_eqsolution cea_solve_eq(
//     const cea_equilibrium_type problem_type,
//     const cea_real state1,
//     const cea_real state2,
//     const cea_int nreactants,
//     const cea_string creactants[],
//     const cea_array amounts,
//     const cea_int nproducts,
//     const cea_string cproducts[],
//     const cea_int ninsert,
//     const cea_string insert[],
//     const char set_trace,
//     const cea_real trace,
//     const char transport,
//     const char ions
// );

//----------------------------------------------------------------------
// Equilibrium Solution API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_eqsolution_create(
    cea_eqsolution *soln,
    const cea_eqsolver solver
);
cea_err cea_eqsolution_destroy(
    cea_eqsolution *soln
);

// Property Queries
cea_err cea_eqsolution_get_property(
    const cea_eqsolution soln,
    const cea_property_type type,
    cea_real *value
);

cea_err cea_eqsolution_get_weights(
    const cea_eqsolution soln,
    const cea_int np,
    cea_real weights[],
    const char log
);

cea_err cea_eqsolution_set_T(
    const cea_eqsolution soln,
    const cea_real T
);

cea_err cea_eqsolution_set_nj(
    const cea_eqsolution soln,
    const cea_eqsolver solver,
    const cea_int np,
    const cea_real nj[]
);

cea_err cea_eqsolution_get_species_amounts(
    const cea_eqsolution soln,
    const cea_int np,
    cea_real amounts[],
    const char mass
);

cea_err cea_eqsolution_get_moles(
    const cea_eqsolution soln,
    cea_real *moles
);

cea_err cea_eqsolution_get_converged(
    const cea_eqsolution soln,
    int *converged
);

//----------------------------------------------------------------------
// Equilibrium Partials API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_eqpartials_create(
    cea_eqpartials *eqpartials,
    const cea_eqsolver solver
);

cea_err cea_eqpartials_destroy(
    cea_eqpartials *eqpartials
);

//----------------------------------------------------------------------
// Rocket Solver API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_rocket_solver_create(
    cea_rocket_solver *solver,
    const cea_mixture products
);

cea_err cea_rocket_solver_create_with_reactants(
    cea_rocket_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants
);

cea_err cea_rocket_solver_create_with_options(
    cea_rocket_solver *solver,
    const cea_mixture products,
    const cea_solver_opts options
);

cea_err cea_rocket_solver_destroy(
    cea_rocket_solver *solver
);

// Get
cea_err cea_rocket_solver_get_size(
    const cea_rocket_solver solver,
    const cea_equilibrium_size eq_variable,
    cea_int *value
);

// Solve
cea_err cea_rocket_solver_solve_iac(
    const cea_rocket_solver solver,
    cea_rocket_solution soln,
    const cea_array weights,
    const cea_real pc,
    const cea_array pi_p,
    const cea_int n_pi_p,
    const cea_array subar,
    const cea_int nsubar,
    const cea_array supar,
    const cea_int nsupar,
    const cea_int n_frz,
    const cea_real hc_or_tc,
    const char use_hc,
    const cea_real tc_est,
    const char use_tc_est
);

cea_err cea_rocket_solver_solve_fac(
    const cea_rocket_solver solver,
    cea_rocket_solution soln,
    const cea_array weights,
    const cea_real pc,
    const cea_array pi_p,
    const cea_int n_pi_p,
    const cea_array subar,
    const cea_int nsubar,
    const cea_array supar,
    const cea_int nsupar,
    const cea_int n_frz,
    const cea_real hc_or_tc,
    const char use_hc,
    const cea_real mdot_or_acat,
    const char use_mdot,
    const cea_real tc_est,
    const char use_tc_est
);

//----------------------------------------------------------------------
// Rocket Solution API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_rocket_solution_create(
    cea_rocket_solution *soln,
    const cea_rocket_solver solver
);

cea_err cea_rocket_solution_destroy(
    cea_rocket_solution *soln
);

// Property Queries
cea_err cea_rocket_solution_get_size(
    const cea_rocket_solution soln,
    cea_int *num_pts
);

cea_err cea_rocket_solution_get_property(
    const cea_rocket_solution soln,
    const cea_rocket_property_type type,
    const cea_int len,
    cea_real *value
);

cea_err cea_rocket_solution_get_weights(
    const cea_rocket_solution soln,
    const cea_int np,
    const cea_int station,
    cea_real weights[],
    const char log
);

cea_err cea_rocket_solution_get_species_amounts(
    const cea_rocket_solution soln,
    const cea_int np,
    const cea_int station,
    cea_real amounts[],
    const char mass
);

cea_err cea_rocket_solution_get_moles(
    const cea_rocket_solution soln,
    cea_real *moles
);

cea_err cea_rocket_solution_get_converged(
    const cea_rocket_solution soln,
    int *converged
);

// cea_err cea_rocket_solution_get_eq_solutions(
//     const cea_rocket_solution soln,
//     cea_int npts,
//     cea_eqsolution *eq_solns[]
// );

// cea_err cea_rocket_solution_destroy_eq_solutions(
//     const cea_rocket_solution soln,
//     cea_int npts,
//     cea_eqsolution *eq_solns[]
// );

//----------------------------------------------------------------------
// Shock Solver API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_shock_solver_create(
    cea_shock_solver *solver,
    const cea_mixture products
);

cea_err cea_shock_solver_create_with_reactants(
    cea_shock_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants
);

cea_err cea_shock_solver_create_with_options(
    cea_shock_solver *solver,
    const cea_mixture products,
    const cea_solver_opts options
);

cea_err cea_shock_solver_destroy(
    cea_shock_solver *solver
);

// Get
cea_err cea_shock_solver_get_size(
    const cea_shock_solver solver,
    const cea_equilibrium_size eq_variable,
    cea_int *value
);

// Solve
cea_err cea_shock_solver_solve(
    const cea_shock_solver solver,
    cea_shock_solution soln,
    const cea_array weights,
    const cea_real T0,
    const cea_real p0,
    const cea_real mach1_or_u1,
    const char use_mach,
    const char refl,
    const char incd_froz,
    const char refl_froz
);

//----------------------------------------------------------------------
// Shock Solution API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_shock_solution_create(
    cea_shock_solution *soln,
    cea_int num_pts
);
cea_err cea_shock_solution_destroy(
    cea_shock_solution *soln
);

// Property Queries
cea_err cea_shock_solution_get_property(
    const cea_shock_solution soln,
    const cea_shock_property_type type,
    const cea_int len,
    cea_real *value
);

cea_err cea_shock_solution_get_scalar_property(
    const cea_shock_solution soln,
    const cea_shock_property_type type,
    cea_real *value
);

cea_err cea_shock_solution_get_weights(
    const cea_shock_solution soln,
    const cea_int np,
    const cea_int station,
    cea_real weights[],
    const char log
);

cea_err cea_shock_solution_get_species_amounts(
    const cea_shock_solution soln,
    const cea_int np,
    const cea_int station,
    cea_real amounts[],
    const char mass
);

cea_err cea_shock_solution_get_moles(
    const cea_shock_solution soln,
    cea_real *moles
);

cea_err cea_shock_solution_get_converged(
    const cea_shock_solution soln,
    int *converged
);

//----------------------------------------------------------------------
// Detonation Solver API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_detonation_solver_create(
    cea_detonation_solver *solver,
    const cea_mixture products
);

cea_err cea_detonation_solver_create_with_reactants(
    cea_detonation_solver *solver,
    const cea_mixture products,
    const cea_mixture reactants
);

cea_err cea_detonation_solver_create_with_options(
    cea_detonation_solver *solver,
    const cea_mixture products,
    const cea_solver_opts options
);

cea_err cea_detonation_solver_destroy(
    cea_detonation_solver *solver
);

// Get
cea_err cea_detonation_solver_get_size(
    const cea_detonation_solver solver,
    const cea_equilibrium_size eq_variable,
    cea_int *value
);

// Solve
cea_err cea_detonation_solver_solve(
    const cea_detonation_solver solver,
    cea_detonation_solution soln,
    const cea_array weights,
    const cea_real T1,
    const cea_real p1,
    const char frozen
);


//----------------------------------------------------------------------
// Detonation Solution API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_detonation_solution_create(
    cea_detonation_solution *soln
);

cea_err cea_detonation_solution_destroy(
    cea_detonation_solution *soln
);

// Property Queries
cea_err cea_detonation_solution_get_property(
    const cea_detonation_solution soln,
    const cea_detonation_property_type type,
    const cea_int len,
    cea_real *value
);

cea_err cea_detonation_solution_get_weights(
    const cea_detonation_solution soln,
    const cea_int np,
    cea_real weights[],
    const char log
);

cea_err cea_detonation_solution_get_species_amounts(
    const cea_detonation_solution soln,
    const cea_int np,
    cea_real amounts[],
    const char mass
);

cea_err cea_detonation_solution_get_moles(
    const cea_detonation_solution soln,
    cea_real *moles
);

cea_err cea_detonation_solution_get_converged(
    const cea_detonation_solution soln,
    int *converged
);

#ifdef __cplusplus
}
#endif
