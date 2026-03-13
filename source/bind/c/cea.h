#include "cea_enum.h"
#include "stdbool.h"

#ifdef __cplusplus
extern "C" {
#endif

// Enumerations
typedef enum { CEA_LOG_LEVEL_ENUM } cea_log_level;
typedef enum { CEA_EQUILIBRIUM_TYPE_ENUM } cea_equilibrium_type;
typedef enum { CEA_DERIVATIVE_METHOD_ENUM } cea_derivative_method;
typedef enum { CEA_EQDERIV_SCALAR_ENUM } cea_eqderiv_scalar;
typedef enum { CEA_EQDERIV_ARRAY_ENUM } cea_eqderiv_array;
typedef enum { CEA_EQDERIV_MATRIX_ENUM } cea_eqderiv_matrix;
typedef enum { CEA_EQUILIBRIUM_SIZE_ENUM } cea_equilibrium_size;
typedef enum { CEA_PROPERTY_TYPE_ENUM } cea_property_type;
typedef enum { CEA_ROCKET_PROPERTY_TYPE_ENUM } cea_rocket_property_type;
typedef enum { CEA_SHOCK_PROPERTY_TYPE_ENUM } cea_shock_property_type;
typedef enum { CEA_DETONATION_PROPERTY_TYPE_ENUM } cea_detonation_property_type;
typedef enum { CEA_ERROR_CODE_ENUM } cea_error_code;

typedef int cea_int;
typedef double cea_real;
typedef double *cea_array;

// Struct-related constants
#define CEA_SPECIES_NAME_LEN 15
#define CEA_ELEMENT_NAME_LEN 2

#ifdef ISO_Fortran_binding_h
#include <ISO_Fortran_binding.h>

// somehow the intel fortran binding header defines CFI_CDESC_T incorrectly
// so we need to undef it and define it ourselves
#undef CFI_CDESC_T
#define CFI_CDESC_T(r)                                                         \
  struct {                                                                     \
    void *base_addr;                                                           \
    size_t elem_len;                                                           \
    size_t reserved1;                                                          \
    size_t flags;                                                              \
    size_t rank;                                                               \
    size_t reserved2;                                                          \
    CFI_dim_t dim[r];                                                          \
  }

struct cea_formula_t {
  CFI_CDESC_T(1) elements;
  CFI_CDESC_T(1) coefficients;
};

struct cea_thermo_fit_t {
  cea_real a1, a2, a3, a4, a5, a6, a7, b1, b2;
};

struct cea_species_thermo_t {
  char name[CEA_SPECIES_NAME_LEN];
  CFI_CDESC_T(0) formula; // type(Formula), allocatable :: formula
  cea_int i_phase;
  cea_int num_intervals;
  cea_real molecular_weight;
  CFI_CDESC_T(2) T_fit;
  CFI_CDESC_T(1) fits;
  cea_real enthalpy_ref;
  cea_real T_ref;
};

struct cea_mixture_t {
  cea_int num_species;
  cea_int num_elements;
  cea_int num_gas;
  cea_int num_condensed;

  CFI_CDESC_T(1) species;      // type(SpeciesThermo), allocatable :: species(:)
  CFI_CDESC_T(1) is_condensed; // logical, allocatable :: is_condensed(:)
  CFI_CDESC_T(2) stoich_matrix; // real(dp), allocatable :: stoich_matrix(:,:)

  CFI_CDESC_T(1)
  species_names; // character(snl), allocatable :: species_names(:)
  CFI_CDESC_T(1)
  element_names; // character(enl), allocatable :: element_names(:)
  bool ions;
};

struct cea_mixture_thermo_t {
  cea_int num_species;
  CFI_CDESC_T(1) cp;
  CFI_CDESC_T(1) cv;
  CFI_CDESC_T(1) enthalpy;
  CFI_CDESC_T(1) entropy;
  CFI_CDESC_T(1) energy;
};

struct cea_eqconstraints_t {
  char type[2];
  cea_real state1;
  cea_real state2;
  CFI_CDESC_T(1) b0;
};

struct cea_eqsolution_t {
  // Inputs
  CFI_CDESC_T(1) w0;

  // Mixture data
  cea_real T;
  CFI_CDESC_T(1) nj;
  cea_real n;
  struct cea_mixture_thermo_t thermo;
  CFI_CDESC_T(1) ln_nj;

  // Solution update variables
  CFI_CDESC_T(1) pi;
  CFI_CDESC_T(1) pi_prev;
  cea_real pi_e;
  cea_real dpi_e;
  CFI_CDESC_T(1) dln_nj;
  CFI_CDESC_T(1) dnj_c;
  cea_real dln_n;
  cea_real dln_T;

  // Algorithm workspace
  CFI_CDESC_T(2) G;
  struct cea_eqconstraints_t constraints;
  CFI_CDESC_T(1) is_active;
  CFI_CDESC_T(1) active_rank;
  cea_int j_liq;
  cea_int j_sol;
  cea_int j_switch;
  cea_int last_cond_idx;
  cea_real T_seed;
  cea_real n_seed;
  CFI_CDESC_T(1) nj_seed;
  CFI_CDESC_T(1) ln_nj_seed;
  CFI_CDESC_T(1) is_active_seed;
  CFI_CDESC_T(1) active_rank_seed;
  cea_int j_liq_seed;
  cea_int j_sol_seed;
  cea_int j_switch_seed;
  cea_int last_cond_idx_seed;

  // Convenience variables
  int32_t gas_converged;
  int32_t condensed_converged;
  int32_t moles_converged;
  int32_t element_converged;
  int32_t temperature_converged;
  int32_t entropy_converged;
  int32_t pi_converged;
  int32_t ions_converged;
  int32_t converged;
  cea_int times_converged;

  // Legacy-style transport component basis cached after convergence.
  cea_int transport_basis_rows;
  CFI_CDESC_T(1) transport_component_idx;
  CFI_CDESC_T(2) transport_basis_matrix;

  // Other solution variables
  CFI_CDESC_T(1) mole_fractions;
  CFI_CDESC_T(1) mass_fractions;

  // Mixture properties
  cea_real density;
  cea_real pressure;
  cea_real volume;
  cea_real M;
  cea_real MW;
  cea_real enthalpy;
  cea_real energy;
  cea_real gibbs_energy;
  cea_real entropy;

  // Transport properties
  cea_real gamma_s;
  cea_real viscosity;
  cea_real cp_fr;
  cea_real cp_eq;
  cea_real cv_fr;
  cea_real cv_eq;
  cea_real conductivity_fr;
  cea_real conductivity_eq;
  cea_real Pr_fr;
  cea_real Pr_eq;
};

struct cea_eqpartials_t {
  CFI_CDESC_T(1) dpi_dlnT;
  CFI_CDESC_T(1) dnc_dlnT;
  cea_real dn_dlnT;
  cea_real dlnV_dlnT;
  CFI_CDESC_T(1) dpi_dlnP;
  CFI_CDESC_T(1) dnc_dlnP;
  cea_real dn_dlnP;
  cea_real dlnV_dlnP;
  cea_real cp_eq;
  cea_real gamma_s;
};

struct cea_eqderivatives_t {
  cea_int m;
  cea_int n;

  // Solver workspace
  CFI_CDESC_T(1) R;
  CFI_CDESC_T(2) J;
  CFI_CDESC_T(2) Rx;
  CFI_CDESC_T(2) dudx;
  CFI_CDESC_T(2) delta_check;

  // Final unpacked derivatives
  cea_real dT_dstate1;
  cea_real dT_dstate2;
  CFI_CDESC_T(1) dT_dw0;

  cea_real dn_dstate1;
  cea_real dn_dstate2;
  CFI_CDESC_T(1) dn_dw0;

  CFI_CDESC_T(1) dnj_dstate1;
  CFI_CDESC_T(1) dnj_dstate2;
  CFI_CDESC_T(2) dnj_dw0;

  cea_real dH_dstate1;
  cea_real dH_dstate2;
  CFI_CDESC_T(1) dH_dw0;

  cea_real dU_dstate1;
  cea_real dU_dstate2;
  CFI_CDESC_T(1) dU_dw0;

  cea_real dG_dstate1;
  cea_real dG_dstate2;
  CFI_CDESC_T(1) dG_dw0;

  cea_real dS_dstate1;
  cea_real dS_dstate2;
  CFI_CDESC_T(1) dS_dw0;

  cea_real dCp_fr_dstate1;
  cea_real dCp_fr_dstate2;
  CFI_CDESC_T(1) dCp_fr_dw0;

  // Finite-difference derivatives (for verification)
  cea_real dT_dstate1_fd;
  cea_real dT_dstate2_fd;
  CFI_CDESC_T(1) dT_dw0_fd;

  cea_real dn_dstate1_fd;
  cea_real dn_dstate2_fd;
  CFI_CDESC_T(1) dn_dw0_fd;

  CFI_CDESC_T(1) dnj_dstate1_fd;
  CFI_CDESC_T(1) dnj_dstate2_fd;
  CFI_CDESC_T(2) dnj_dw0_fd;

  cea_real dH_dstate1_fd;
  cea_real dH_dstate2_fd;
  CFI_CDESC_T(1) dH_dw0_fd;

  cea_real dU_dstate1_fd;
  cea_real dU_dstate2_fd;
  CFI_CDESC_T(1) dU_dw0_fd;

  cea_real dG_dstate1_fd;
  cea_real dG_dstate2_fd;
  CFI_CDESC_T(1) dG_dw0_fd;

  cea_real dS_dstate1_fd;
  cea_real dS_dstate2_fd;
  CFI_CDESC_T(1) dS_dw0_fd;

  cea_real dCp_fr_dstate1_fd;
  cea_real dCp_fr_dstate2_fd;
  CFI_CDESC_T(1) dCp_fr_dw0_fd;
};

#endif

// Type Shorthand, Incomplete types
typedef struct cea_reactant_t *cea_reactant;
typedef struct cea_mixture_t *cea_mixture;
typedef struct cea_eqsolver_t *cea_eqsolver;
typedef struct cea_eqsolution_t *cea_eqsolution;
typedef struct cea_eqpartials_t *cea_eqpartials;
typedef struct cea_eqderivatives_t *cea_eqderivatives;
typedef struct cea_rocket_solver_t *cea_rocket_solver;
typedef struct cea_rocket_solution_t *cea_rocket_solution;
typedef struct cea_shock_solver_t *cea_shock_solver;
typedef struct cea_shock_solution_t *cea_shock_solution;
typedef struct cea_detonation_solver_t *cea_detonation_solver;
typedef struct cea_detonation_solution_t *cea_detonation_solution;
typedef cea_error_code cea_err;
typedef const char *cea_string;

typedef struct {
  cea_string name;
  cea_int num_elements;
  const cea_string *elements;
  const cea_real *coefficients;
  bool has_molecular_weight;
  cea_real molecular_weight;
  bool has_enthalpy;
  cea_real enthalpy;
  cea_string enthalpy_units;
  bool has_temperature;
  cea_real temperature;
  cea_string temperature_units;
} cea_reactant_input;

// Struct types for optional arguments
typedef struct {
  cea_real trace;
  bool ions;
  bool transport;
  cea_mixture reactants;
  cea_int ninsert;
  const cea_string *insert;
  bool smooth_truncation; // enable smooth logistic truncation instead of hard
                          // cutoff (default false)
  cea_real truncation_width; // gate width in log-space; <= 0 means use solver
                             // default (0.25)
} cea_solver_opts;

// Initialize optional arguments
cea_err cea_solver_opts_init(cea_solver_opts *opts);

// Input string ownership
// All input strings/arrays (species, reactants, omit, insert, reactant_input
// fields) are copied by the C/Fortran layer during the call. Callers retain
// ownership and may free or reuse their buffers after the call returns.

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
cea_err cea_mixture_create(cea_mixture *mix, const cea_int nspecies,
                           const cea_string species[]);

cea_err cea_mixture_create_w_ions(cea_mixture *mix, const cea_int nspecies,
                                  const cea_string species[]);

cea_err cea_mixture_create_from_reactants(cea_mixture *mix,
                                          const cea_int nreactants,
                                          const cea_string reactants[],
                                          const cea_int nomit,
                                          const cea_string omit[]);

cea_err cea_mixture_create_from_reactants_w_ions(cea_mixture *mix,
                                                 const cea_int nreactants,
                                                 const cea_string reactants[],
                                                 const cea_int nomit,
                                                 const cea_string omit[]);

cea_err
cea_mixture_create_from_input_reactants(cea_mixture *mix,
                                        const cea_int nreactants,
                                        const cea_reactant_input reactants[]);

cea_err cea_mixture_create_from_input_reactants_w_ions(
    cea_mixture *mix, const cea_int nreactants,
    const cea_reactant_input reactants[]);

cea_err cea_mixture_create_products_from_input_reactants(
    cea_mixture *mix, const cea_int nreactants,
    const cea_reactant_input reactants[], const cea_int nomit,
    const cea_string omit[]);

cea_err cea_mixture_create_products_from_input_reactants_w_ions(
    cea_mixture *mix, const cea_int nreactants,
    const cea_reactant_input reactants[], const cea_int nomit,
    const cea_string omit[]);

cea_err cea_mixture_destroy(cea_mixture *mix);

// Get values
cea_err cea_mixture_get_num_species(const cea_mixture mix,
                                    cea_int *num_species);

// Deprecated: returns a heap-allocated string that must be freed with
// cea_string_free.
cea_err cea_mixture_get_species_name(const cea_mixture *mix,
                                     const cea_int i_species,
                                     cea_string *species);

// Deprecated: returns heap-allocated strings that must be freed with
// cea_string_array_free.
cea_err cea_mixture_get_species_names(const cea_mixture *mix,
                                      const cea_int nspecies,
                                      cea_string *species[]);

// Buffer-based species name retrieval
cea_err cea_mixture_get_species_name_buf(const cea_mixture *mix,
                                         const cea_int i_species, char *species,
                                         const cea_int buf_len);

cea_err cea_mixture_get_species_names_buf(const cea_mixture *mix,
                                          const cea_int nspecies, char *species,
                                          const cea_int stride);

// Free functions for deprecated allocating getters
cea_err cea_string_free(cea_string species);
cea_err cea_string_array_free(cea_string species[], const cea_int nspecies);

// Unit Conversion
cea_err cea_mixture_moles_to_weights(const cea_mixture mix, const cea_int len,
                                     const cea_real moles[],
                                     cea_real weights[]);

cea_err cea_mixture_weights_to_moles(const cea_mixture mix, const cea_int len,
                                     const cea_real weights[],
                                     cea_real moles[]);

cea_err cea_mixture_per_mole_to_per_weight(const cea_mixture mix,
                                           const cea_int len,
                                           const cea_real per_mole[],
                                           cea_real per_weight[]);

cea_err cea_mixture_per_weight_to_per_mole(const cea_mixture mix,
                                           const cea_int len,
                                           const cea_real per_weight[],
                                           cea_real per_mole[]);

// Fuel/Oxidant Mixture Conversion
cea_err cea_mixture_chem_eq_ratio_to_of_ratio(const cea_mixture mix,
                                              const cea_int len,
                                              const cea_real oxidant_weights[],
                                              const cea_real fuel_weights[],
                                              const cea_real chem_eq_ratio,
                                              cea_real *of_ratio);

cea_err cea_mixture_weight_eq_ratio_to_of_ratio(
    const cea_mixture mix, const cea_int len, const cea_real oxidant_weights[],
    const cea_real fuel_weights[], const cea_real weight_eq_ratio,
    cea_real *of_ratio);

cea_err cea_mixture_of_ratio_to_weights(const cea_mixture mix,
                                        const cea_int len,
                                        const cea_real oxidant_weights[],
                                        const cea_real fuel_weights[],
                                        const cea_real of_ratio,
                                        cea_real reactant_weights[]);

// Property Calculation
cea_err cea_mixture_calc_property(const cea_mixture mix,
                                  const cea_property_type type,
                                  const cea_int len_weights,
                                  const cea_real weights[],
                                  const cea_real temperature, cea_real *value);

cea_err cea_mixture_calc_property_multitemp(const cea_mixture mix,
                                            const cea_property_type type,
                                            const cea_int len_weights,
                                            const cea_real weights[],
                                            const cea_int len_temperatures,
                                            const cea_real temperatures[],
                                            cea_real *value);

cea_err cea_mixture_calc_property_tp(const cea_mixture mix,
                                     const cea_property_type type,
                                     const cea_int len_weights,
                                     const cea_real weights[],
                                     const cea_real temperature,
                                     const cea_real pressure, cea_real *value);

cea_err cea_mixture_calc_property_tp_multitemp(
    const cea_mixture mix, const cea_property_type type,
    const cea_int len_weights, const cea_real weights[],
    const cea_int len_temperatures, const cea_real temperatures[],
    const cea_real pressure, cea_real *value);

//----------------------------------------------------------------------
// Equilibrium Solver API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_eqsolver_create(cea_eqsolver *solver, const cea_mixture products);

cea_err cea_eqsolver_create_with_reactants(cea_eqsolver *solver,
                                           const cea_mixture products,
                                           const cea_mixture reactants);

cea_err cea_eqsolver_create_with_options(cea_eqsolver *solver,
                                         const cea_mixture products,
                                         const cea_solver_opts options);

cea_err cea_eqsolver_destroy(cea_eqsolver *solver);

// Solve
cea_err cea_eqsolver_solve(const cea_eqsolver solver,
                           const cea_equilibrium_type type,
                           const cea_real state1, const cea_real state2,
                           const cea_real amounts[], cea_eqsolution soln);

cea_err cea_eqsolver_solve_with_partials(
    const cea_eqsolver solver, const cea_equilibrium_type type,
    const cea_real state1, const cea_real state2, const cea_real amounts[],
    cea_eqsolution soln, cea_eqpartials eqpartials);

// Querry functions
cea_err cea_eqsolver_get_size(const cea_eqsolver solver,
                              const cea_equilibrium_size eq_variable,
                              cea_int *value);

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
//     const bool set_trace,
//     const cea_real trace,
//     const bool transport,
//     const bool ions
// );

//----------------------------------------------------------------------
// Equilibrium Solution API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_eqsolution_create(cea_eqsolution *soln, const cea_eqsolver solver);
cea_err cea_eqsolution_destroy(cea_eqsolution *soln);

// Property Queries
cea_err cea_eqsolution_get_property(const cea_eqsolution soln,
                                    const cea_property_type type,
                                    cea_real *value);

cea_err cea_eqsolution_get_weights(const cea_eqsolution soln, const cea_int np,
                                   cea_real weights[], const bool log);

cea_err cea_eqsolution_set_T(const cea_eqsolution soln, const cea_real T);

cea_err cea_eqsolution_set_nj(const cea_eqsolution soln,
                              const cea_eqsolver solver, const cea_int np,
                              const cea_real nj[]);

cea_err cea_eqsolution_get_species_amounts(const cea_eqsolution soln,
                                           const cea_int np, cea_real amounts[],
                                           const bool mass);

cea_err cea_eqsolution_get_moles(const cea_eqsolution soln, cea_real *moles);

cea_err cea_eqsolution_get_converged(const cea_eqsolution soln, int *converged);

//----------------------------------------------------------------------
// Equilibrium Partials API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_eqpartials_create(cea_eqpartials *eqpartials,
                              const cea_eqsolver solver);

cea_err cea_eqpartials_destroy(cea_eqpartials *eqpartials);

//----------------------------------------------------------------------
// Equilibrium Derivatives API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_eqderivatives_create(cea_eqderivatives *derivs,
                                 const cea_eqsolver solver,
                                 const cea_eqsolution soln);

cea_err cea_eqderivatives_destroy(cea_eqderivatives *derivs);

// Compute
cea_err cea_eqderivatives_compute_derivatives(const cea_eqderivatives derivs,
                                              const cea_eqsolver solver,
                                              const cea_eqsolution soln,
                                              const bool check_closure_defect);

cea_err cea_eqderivatives_compute_fd(const cea_eqderivatives derivs,
                                     const cea_eqsolver solver,
                                     const cea_eqsolution soln,
                                     const cea_real h, const bool verbose,
                                     const bool central);

// Getters
cea_err cea_eqderivatives_get_scalar(const cea_eqderivatives derivs,
                                     const cea_eqderiv_scalar which,
                                     const cea_derivative_method method,
                                     cea_real *value);

cea_err cea_eqderivatives_get_array(const cea_eqderivatives derivs,
                                    const cea_eqsolver solver,
                                    const cea_eqsolution soln,
                                    const cea_eqderiv_array which,
                                    const cea_derivative_method method,
                                    const cea_int len, cea_real out[]);

// Matrix is returned row-major (C-order) in out[(row*cols) + col].
cea_err cea_eqderivatives_get_matrix(const cea_eqderivatives derivs,
                                     const cea_eqsolver solver,
                                     const cea_eqsolution soln,
                                     const cea_eqderiv_matrix which,
                                     const cea_derivative_method method,
                                     const cea_int rows, const cea_int cols,
                                     cea_real out[]);

//----------------------------------------------------------------------
// Rocket Solver API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_rocket_solver_create(cea_rocket_solver *solver,
                                 const cea_mixture products);

cea_err cea_rocket_solver_create_with_reactants(cea_rocket_solver *solver,
                                                const cea_mixture products,
                                                const cea_mixture reactants);

cea_err cea_rocket_solver_create_with_options(cea_rocket_solver *solver,
                                              const cea_mixture products,
                                              const cea_solver_opts options);

cea_err cea_rocket_solver_destroy(cea_rocket_solver *solver);

// Get
cea_err cea_rocket_solver_get_size(const cea_rocket_solver solver,
                                   const cea_equilibrium_size eq_variable,
                                   cea_int *value);

// Solve
cea_err cea_rocket_solver_solve_iac(
    const cea_rocket_solver solver, cea_rocket_solution soln,
    const cea_real weights[], const cea_real pc,
    // Optional: ignored when n_pi_p == 0 (pi_p may be NULL)
    const cea_real pi_p[], const cea_int n_pi_p, const cea_real subar[],
    const cea_int nsubar, const cea_real supar[], const cea_int nsupar,
    const cea_int n_frz, const cea_real hc_or_tc, const bool use_hc,
    const cea_real tc_est, const bool use_tc_est);

cea_err cea_rocket_solver_solve_fac(
    const cea_rocket_solver solver, cea_rocket_solution soln,
    const cea_real weights[], const cea_real pc,
    // Optional: ignored when n_pi_p == 0 (pi_p may be NULL)
    const cea_real pi_p[], const cea_int n_pi_p, const cea_real subar[],
    const cea_int nsubar, const cea_real supar[], const cea_int nsupar,
    const cea_int n_frz, const cea_real hc_or_tc, const bool use_hc,
    const cea_real mdot_or_acat, const bool use_mdot, const cea_real tc_est,
    const bool use_tc_est);

//----------------------------------------------------------------------
// Rocket Solution API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_rocket_solution_create(cea_rocket_solution *soln,
                                   const cea_rocket_solver solver);

cea_err cea_rocket_solution_destroy(cea_rocket_solution *soln);

// Property Queries
cea_err cea_rocket_solution_get_size(const cea_rocket_solution soln,
                                     cea_int *num_pts);

cea_err cea_rocket_solution_get_property(const cea_rocket_solution soln,
                                         const cea_rocket_property_type type,
                                         const cea_int len, cea_real *value);

cea_err cea_rocket_solution_get_weights(const cea_rocket_solution soln,
                                        const cea_int np, const cea_int station,
                                        cea_real weights[], const bool log);

cea_err cea_rocket_solution_get_species_amounts(const cea_rocket_solution soln,
                                                const cea_int np,
                                                const cea_int station,
                                                cea_real amounts[],
                                                const bool mass);

cea_err cea_rocket_solution_get_moles(const cea_rocket_solution soln,
                                      cea_real *moles);

cea_err cea_rocket_solution_get_converged(const cea_rocket_solution soln,
                                          int *converged);

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
cea_err cea_shock_solver_create(cea_shock_solver *solver,
                                const cea_mixture products);

cea_err cea_shock_solver_create_with_reactants(cea_shock_solver *solver,
                                               const cea_mixture products,
                                               const cea_mixture reactants);

cea_err cea_shock_solver_create_with_options(cea_shock_solver *solver,
                                             const cea_mixture products,
                                             const cea_solver_opts options);

cea_err cea_shock_solver_destroy(cea_shock_solver *solver);

// Get
cea_err cea_shock_solver_get_size(const cea_shock_solver solver,
                                  const cea_equilibrium_size eq_variable,
                                  cea_int *value);

// Solve
cea_err cea_shock_solver_solve(const cea_shock_solver solver,
                               cea_shock_solution soln,
                               const cea_real weights[], const cea_real T0,
                               const cea_real p0, const cea_real mach1_or_u1,
                               const bool use_mach, const bool refl,
                               const bool incd_froz, const bool refl_froz);

//----------------------------------------------------------------------
// Shock Solution API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_shock_solution_create(cea_shock_solution *soln, cea_int num_pts);
cea_err cea_shock_solution_destroy(cea_shock_solution *soln);

// Property Queries
cea_err cea_shock_solution_get_property(const cea_shock_solution soln,
                                        const cea_shock_property_type type,
                                        const cea_int len, cea_real *value);

cea_err
cea_shock_solution_get_scalar_property(const cea_shock_solution soln,
                                       const cea_shock_property_type type,
                                       cea_real *value);

cea_err cea_shock_solution_get_weights(const cea_shock_solution soln,
                                       const cea_int np, const cea_int station,
                                       cea_real weights[], const bool log);

cea_err cea_shock_solution_get_species_amounts(const cea_shock_solution soln,
                                               const cea_int np,
                                               const cea_int station,
                                               cea_real amounts[],
                                               const bool mass);

cea_err cea_shock_solution_get_moles(const cea_shock_solution soln,
                                     cea_real *moles);

cea_err cea_shock_solution_get_converged(const cea_shock_solution soln,
                                         int *converged);

//----------------------------------------------------------------------
// Detonation Solver API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_detonation_solver_create(cea_detonation_solver *solver,
                                     const cea_mixture products);

cea_err
cea_detonation_solver_create_with_reactants(cea_detonation_solver *solver,
                                            const cea_mixture products,
                                            const cea_mixture reactants);

cea_err
cea_detonation_solver_create_with_options(cea_detonation_solver *solver,
                                          const cea_mixture products,
                                          const cea_solver_opts options);

cea_err cea_detonation_solver_destroy(cea_detonation_solver *solver);

// Get
cea_err cea_detonation_solver_get_size(const cea_detonation_solver solver,
                                       const cea_equilibrium_size eq_variable,
                                       cea_int *value);

// Solve
cea_err cea_detonation_solver_solve(const cea_detonation_solver solver,
                                    cea_detonation_solution soln,
                                    const cea_real weights[], const cea_real T1,
                                    const cea_real p1, const bool frozen);

//----------------------------------------------------------------------
// Detonation Solution API
//----------------------------------------------------------------------

// Create/Destroy
cea_err cea_detonation_solution_create(cea_detonation_solution *soln);

cea_err cea_detonation_solution_destroy(cea_detonation_solution *soln);

// Property Queries
cea_err
cea_detonation_solution_get_property(const cea_detonation_solution soln,
                                     const cea_detonation_property_type type,
                                     const cea_int len, cea_real *value);

cea_err cea_detonation_solution_get_weights(const cea_detonation_solution soln,
                                            const cea_int np,
                                            cea_real weights[], const bool log);

cea_err cea_detonation_solution_get_species_amounts(
    const cea_detonation_solution soln, const cea_int np, cea_real amounts[],
    const bool mass);

cea_err cea_detonation_solution_get_moles(const cea_detonation_solution soln,
                                          cea_real *moles);

cea_err
cea_detonation_solution_get_converged(const cea_detonation_solution soln,
                                      int *converged);

#ifdef __cplusplus
}
#endif
