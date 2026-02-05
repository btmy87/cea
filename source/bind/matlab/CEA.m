classdef CEA
    properties
        alias = "cea";
        warnFailure = true;
        errFailure = false;
    end
    
    methods
        function obj = CEA(opts)
            arguments
                opts.alias (1, 1) string = "cea";
                opts.warnFailure(1,1) {mustBeNumericOrLogical} = true
                opts.errFailure(1,1) {mustBeNumericOrLogical} = false
            end
            % Create 
            if ~libisloaded(obj.alias)
                load_cea(alias=opts.alias);
            end

            obj.warnFailure = opts.warnFailure;
            obj.errFailure = obj.errFailure;
        end

        function checkval(obj, value, msg)
            arguments
                obj
                value (1, 1) string
                msg (1, 1) string = ""
            end
            % checks cea error values
            if value ~= "CEA_SUCCESS"
                if ~isempty(msg) && msg ~= ""
                    fprintf("%s\n", msg);
                end
                if obj.errFailure
                    error("CEA call failed: %s", value)
                elseif obj.warnFailure
                    warning("CEA call failed: %s", value)
                end
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
            obj.checkval(err);
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
            obj.checkval(err, "Error in cea_mixture_create");
        end
        function mix = mixture_create_from_reactants(obj, reactants, omitted)
            arguments
                obj
                reactants (1, :) string
                omitted (1, :) string
            end
            nr = length(reactants);
            no = length(omitted);
            mix = libpointer("voidPtr", 0);
            reac = libpointer("stringPtrPtr", cellstr(reactants));
            omit = libpointer("stringPtrPtr", cellstr(omitted));
            err = calllib(obj.alias, "cea_mixture_create_from_reactants", ...
                    mix, nr, reac, no, omit);
            obj.checkval(err, "Error in cea_mixture_create_from_reactants");
        end

        function out = mixture_calc_property_multitemp(obj, mix, type, weights, temperatures)
            outPtr = libpointer("doublePtr", 0);
            err = calllib(obj.alias, "cea_mixture_calc_property_multitemp", ...
                    mix, type, length(weights), weights, ...
                    length(temperatures), temperatures, outPtr);
            out = outPtr.Value;

            obj.checkval(err, "Error in cea_mixture_calc_property_multitemp");
        end
           
        function opts = solver_opts_init(obj)
            opts = libpointer("cea_solver_opts", struct());
            err = calllib(obj.alias, "cea_solver_opts_init", opts);
            obj.checkval(err, "Error in cea_solver_opts_init");
        end

        function solver = eqsolver_create_with_options(obj, prod, opts)
            solver = libpointer("voidPtr", 0);
            err = calllib(obj.alias, "cea_eqsolver_create_with_options", ...
                solver, prod, opts);
            obj.checkval(err, "Error in eqsolver_create_with_options");
        end

        function solver = eqsolver_create_with_reactants(obj, prod, reac)
            solver = libpointer("voidPtr", 0);
            err = calllib(obj.alias, "cea_eqsolver_create_with_reactants", ...
                solver, prod, reac);
            obj.checkval(err, "Error in eqsolver_create_with_reactants");
        end

        function solution = eqsolution_create(obj, solver)
            solution = libpointer("voidPtr", 0);
            err = calllib(obj.alias, "cea_eqsolution_create", ...
                solution, solver);
            obj.checkval(err, "Error in cea_eqsolution_create");
        end

        function partials = eqpartials_create(obj, solver)
            partials = libpointer("voidPtr", 0);
            err = calllib(obj.alias, "cea_eqpartials_create", ...
                partials, solver);
            obj.checkval(err, "Error in cea_eqpartialscreate");
        end

        function weights = mixture_moles_to_weights(obj, reac, n, moles)
            weights = libpointer("doublePtr", zeros(1, n));
            err = calllib(obj.alias, "cea_mixture_moles_to_weights", ...
                reac, n, moles, weights);
            obj.checkval(err, "Error in cea_mixture_moles_to_weights");
        end

        function OFRatio = mixture_chem_eq_ratio_to_of_ratio(...
                obj, reac, n, oxidWeights, fuelWeights, er)
            ofptr = libpointer("doublePtr", 0);
            err = calllib(obj.alias, "cea_mixture_chem_eq_ratio_to_of_ratio", ...
                reac, n, oxidWeights, fuelWeights, er, ofptr);
            obj.checkval(err, "Error in cea_mixture_chem_eq_ratio_to_of_ratio");
            OFRatio = ofptr.Value;
        end

        function weights = mixture_of_ratio_to_weights(...
                obj, reac, n, oxidWeights, fuelWeights, OFRatio)
            weightsPtr = libpointer("doublePtr", zeros(1, n));
            err = calllib(obj.alias, "cea_mixture_of_ratio_to_weights", ...
                reac, n, oxidWeights, fuelWeights, OFRatio, weightsPtr);
            obj.checkval(err, "Error in cea_mixture_of_ratio_to_weights");
            weights = weightsPtr.Value;
        end

        function eqsolver_solve(obj, solver, prob, ...
                state1, state2, weights, soln)
            err = calllib(obj.alias, "cea_eqsolver_solve", ...
                solver, prob, state1, state2, weights, soln);
            obj.checkval(err, "Error in cea_eqsolver_solve");
        end

        function eqsolver_solve_with_partials(obj, solver, prob, ...
                state1, state2, weights, soln, partials)
            err = calllib(obj.alias, "cea_eqsolver_solve_with_partials", ...
                solver, prob, state1, state2, weights, soln, partials);
            obj.checkval(err, "Error in cea_eqsolver_solve_with_partials");
        end

        function out = eqsolution_get_property(obj, soln, prop)
            outPtr = libpointer("doublePtr", 0.0);
            err = calllib(obj.alias, "cea_eqsolution_get_property", ...
                soln, prop, outPtr);
            obj.checkval(err, "Error in cea_eqsolution_get_property");
            out = outPtr.Value; 
        end

        function eqpartials_destroy(obj, partials)
            err = calllib(obj.alias, "cea_eqpartials_destroy", partials);
            obj.checkval(err, "Error in cea_eqpartials_destroy");
        end

        function eqsolution_destroy(obj, soln)
            err = calllib(obj.alias, "cea_eqsolution_destroy", soln);
            obj.checkval(err, "Error in cea_eqsolution_destroy");
        end

        function eqsolver_destroy(obj, solver)
            err = calllib(obj.alias, "cea_eqsolver_destroy", solver);
            obj.checkval(err, "Error in cea_eqsolver_destroy");
        end

        function mixture_destroy(obj, mix)
            err = calllib(obj.alias, "cea_mixture_destroy", mix);
            obj.checkval(err, "Error in cea_mixture_destroy");
        end

    end


end
