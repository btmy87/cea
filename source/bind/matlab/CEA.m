classdef CEA
    properties
        alias = "cea";
    
    end
    
    methods
        function obj = CEA()
            % Create 
            if ~libisloaded(obj.alias)
                load_cea();
            end
        end

        function set_log_level(obj, logLevel)
            arguments
                obj
                logLevel (1, 1) string {mustBeMember(logLevel, ...
                    ["CEA_LOG_CRITICAL", "CEA_LOG_ERROR", ...
                     "CEA_LOG_WARNING", "CEA_LOG_INFO", ...
                     "CEA_LOG_DEBUG", "CEA_LOG_NONE"])} = ...
                     "CEA_LOG_CRITICAL";
            end
            err = calllib(obj.alias, "cea_set_log_level", char(logLevel));
            assert(err == "CEA_SUCCESS");
        end

        function mix = mixture_create(obj, reactants)
            arguments
                obj
                reactants (1, :) string;
            end
            mix = libpointer("voidPtr", 0);
            n = length(reactants);
            % reac = libpointer("cstring", cellstr(reactants));
            reac = libpointer("stringPtrPtr", cellstr(reactants));
            err = calllib(obj.alias, "cea_mixture_create", ...
                mix, n, reac);
            assert(err == "CEA_SUCCESS", "Error in cea_mixture_create");
        end
           
        function opts = solver_opts_init(obj)
            opts = libpointer("cea_solver_opts", struct());
            err = calllib(obj.alias, "cea_solver_opts_init", opts);
            assert(err == "CEA_SUCCESS", "Error in cea_solver_opts_init");
        end

        function solver = eqsolver_create_with_options(obj, prod, opts)
            solver = libpointer("voidPtr", 0);
            err = calllib(obj.alias, "cea_eqsolver_create_with_options", ...
                solver, prod, opts);
            assert(err == "CEA_SUCCESS", "Error in eqsolver_create_with_options");
        end

        function solution = eqsolution_create(obj, solver)
            solution = libpointer("voidPtr", 0);
            err = calllib(obj.alias, "cea_eqsolution_create", ...
                solution, solver);
            assert(err == "CEA_SUCCESS", "Error in cea_eqsolution_create");
        end

        function partials = eqpartials_create(obj, solver)
            partials = libpointer("voidPtr", 0);
            err = calllib(obj.alias, "cea_eqpartials_create", ...
                partials, solver);
            assert(err == "CEA_SUCCESS", "Error in cea_eqpartialscreate");
        end

        function weights = mixture_moles_to_weights(obj, reac, n, moles)
            weights = libpointer("doublePtr", zeros(1, n));
            err = calllib(obj.alias, "cea_mixture_moles_to_weights", ...
                reac, n, moles, weights);
            assert(err == "CEA_SUCCESS", "Error in cea_mixture_moles_to_weights");
        end

        function OFRatio = mixture_chem_eq_ratio_to_of_ratio(...
                obj, reac, n, oxidWeights, fuelWeights, er)
            ofptr = libpointer("doublePtr", 0);
            err = calllib(obj.alias, "cea_mixture_chem_eq_ratio_to_of_ratio", ...
                reac, n, oxidWeights, fuelWeights, er, ofptr);
            assert(err == "CEA_SUCCESS", "Error in cea_mixture_chem_eq_ratio_to_of_ratio");
            OFRatio = ofptr.Value;
        end

        function weights = mixture_of_ratio_to_weights(...
                obj, reac, n, oxidWeights, fuelWeights, OFRatio)
            weightsPtr = libpointer("doublePtr", zeros(1, n));
            err = calllib(obj.alias, "cea_mixture_of_ratio_to_weights", ...
                reac, n, oxidWeights, fuelWeights, OFRatio, weightsPtr);
            assert(err == "CEA_SUCCESS", "Error in cea_mixture_of_ratio_to_weights");
            weights = weightsPtr.Value;
        end

        function eqsolver_solve_with_partials(obj, solver, prob, ...
                state1, state2, weights, soln, partials)
            err = calllib(obj.alias, "cea_eqsolver_solve_with_partials", ...
                solver, prob, state1, state2, weights, soln, partials);
            assert(err == "CEA_SUCCESS", "Error in cea_eqsolver_solve_with_partials");
        end

        function out = eqsolution_get_property(obj, soln, prop)
            outPtr = libpointer("doublePtr", 0.0);
            err = calllib(obj.alias, "cea_eqsolution_get_property", ...
                soln, prop, outPtr);
            assert(err == "CEA_SUCCESS", "Error in cea_eqsolution_get_property");
            out = outPtr.Value; 
        end

        function eqpartials_destroy(obj, partials)
            err = calllib(obj.alias, "cea_eqpartials_destroy", partials);
            assert(err == "CEA_SUCCESS", "Error in cea_eqpartials_destroy");
        end

        function eqsolution_destroy(obj, soln)
            err = calllib(obj.alias, "cea_eqsolution_destroy", soln);
            assert(err == "CEA_SUCCESS", "Error in cea_eqsolution_destroy");
        end

        function eqsolver_destroy(obj, solver)
            err = calllib(obj.alias, "cea_eqsolver_destroy", solver);
            assert(err == "CEA_SUCCESS", "Error in cea_eqsolver_destroy");
        end

        function mixture_destroy(obj, mix)
            err = calllib(obj.alias, "cea_mixture_destroy", mix);
            assert(err == "CEA_SUCCESS", "Error in cea_mixture_destroy");
        end
    end


end