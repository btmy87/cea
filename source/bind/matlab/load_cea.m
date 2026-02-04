function load_cea(opts)
% load_cea loads cea bindc dll
%
% INPUT OPTIONS:
%  ceaDataDir: directory containing thermo.lib and trans.lib
%  ceaInstallDir: CEA installation directory, should contain
%                 - lib\cea_bindc.dll
%                 - include\cea\bindc\cea.h
%                 - include\cea\bindc\cea_enum.h
%                 If built and installed, matches CMAKE_INSTALL_PREFIX
%  alias: alias for loading dll, defaults to "cea"
%  makeThunk: Set to true to generate a reusable thunk file
%  useThunk: Set to true to use an existing thunk file generated with
%            makeThunk.  This is faster than recompiling the interface and
%            can be used without a compiler setup.
%  reload: Set to true to force reload if already loaded.  Defaults to
%          false which will do nothing if CEA is already loaded.
%  verbose: Set to true for additional output about loading process
%  shortname: Set to library name, defaults to "cea_bindc", which should
%             work for most cases.
%  initialize: Defaults to true, set to false to skip CEA initialization.
%              Use false to accomodate non-standard library names or
%              locations.  The user will need to call cea_init_thermo and
%              cea_init_trans manually.
%

arguments
    opts.ceaDataDir (1, 1) string = string(getenv("CEA_DATA_DIR"));
    opts.ceaInstallDir (1, 1) string = "";
    opts.alias (1, 1) string = "cea";
    opts.makeThunk (1, 1) {mustBeNumericOrLogical} = false;
    opts.useThunk (1, 1) {mustBeNumericOrLogical} = false;
    opts.reload (1, 1) {mustBeNumericOrLogical} = false;
    opts.verbose (1, 1) {mustBeNumericOrLogical} = false;
    opts.shortname (1, 1) string = "cea_bindc";
    opts.initialize (1, 1) {mustBeNumericOrLogical} = true;
end

%% Check if already loaded, unload if requested
if libisloaded(opts.alias)
    if opts.reload
        unloadlibrary(opts.alias);
    else
        if opts.verbose
            fprintf("CEA already loaded.  Use reload=true to force reload.\n");
        end
        return;
    end
end

%% Validate Iputs
assert(isfolder(opts.ceaDataDir), ...
    "ceaDataDir must be a valid folder: %s", opts.ceaDataDir)

thermoLib = fullfile(opts.ceaDataDir, "thermo.lib");
transLib = fullfile(opts.ceaDataDir, "trans.lib");
assert(isfile(thermoLib), "Couldn't find thermo.lib: %s", thermoLib);
assert(isfile(transLib), "couldn't find trans.lib: %s", transLib);

if isempty(opts.ceaInstallDir) || opts.ceaInstallDir == ""
    opts.ceaInstallDir = fullfile(opts.ceaDataDir, "..");
end

%% Find files

if ispc
    libname = fullfile(opts.ceaInstallDir, "lib", ...
        opts.shortname + ".dll");
else
    libname = fullfile(opts.ceaInstallDir, "lib", ...
        "lib" + opts.shortname + ".so");
end
assert(isfile(libname), "CEA Library does not exist: %s", libname);

incdir = fullfile(opts.ceaInstallDir, "include", "cea", "bindc");
header1 = fullfile(incdir, "cea.h");
header2 = fullfile(incdir, "cea_enum.h");

thunkfile = opts.shortname + "_thunk_" + computer;
protofile = opts.shortname + "_proto_" + computer;


%% Call load_library
if opts.makeThunk
    % if requested, generate a thunk file   
    if opts.verbose
        fprintf("Loading CEA and generating thunkfile.\n");
        fprintf("  libname: %s\n", libname)
        fprintf("  header1: %s\n", header1);
        fprintf("  header2: %s\n", header2);
        fprintf("  thunkfile: %s\n", thunkfile);
        fprintf("  protofile: %s\n", protofile);
    end
    loadlibrary(libname, header1, ...
        addheader=header2, ...
        alias=opts.alias, ...
        mfilename=protofile, ...
        thunkfilename=thunkfile);

    % clean up unneeded outputs
    if ispc
        delete(thunkfile + ".exp");
        delete(thunkfile + ".lib");
        delete(thunkfile + ".obj");
    end
elseif opts.useThunk
    % load from pre-existing thunk file
    if opts.verbose
        fprintf("Loading CEA from thunkfile.\n");
        fprintf("  libname: %s\n", libname)
        fprintf("  thunkfile: %s\n", thunkfile);
        fprintf("  protofile: %s\n", protofile);
    end
    loadlibrary(libname, str2func(protofile), alias=opts.alias);

else
    % load directly
    if opts.verbose
        fprintf("Loading CEA without thunk file.\n");
        fprintf("  libname: %s\n", libname)
        fprintf("  header1: %s\n", header1);
        fprintf("  header2: %s\n", header2);
    end

    assert(isfile(header1), "CEA header1 does not exist: %s", header1);
    assert(isfile(header1), "CEA header2 does not exist: %s", header2);
    loadlibrary(libname, header1, ...
        addheader=header2, ...
        alias=opts.alias)
end

%% Test Library
v1 = libpointer("int32Ptr", 0);
v2 = libpointer("int32Ptr", 0);
v3 = libpointer("int32Ptr", 0);

out1 = calllib(opts.alias, "cea_version_major", v1);
out2 = calllib(opts.alias, "cea_version_minor", v2);
out3 = calllib(opts.alias, "cea_version_patch", v3);
fprintf("Loaded CEA v%d.%d.%d\n", v1.Value, v2.Value, v3.Value);

assert(out1 == "CEA_SUCCESS", "Error calling cea_version_major");
assert(out2 == "CEA_SUCCESS", "Error calling cea_version_minor");
assert(out3 == "CEA_SUCCESS", "Error calling cea_version_patch");

%% Initialize library
if opts.initialize
    thermolib = fullfile(opts.ceaDataDir, "thermo.lib");
    translib = fullfile(opts.ceaDataDir, "trans.lib");
    if opts.verbose
        fprintf("Initializing CEA thermo: %s\n", thermolib);
    end
    calllib(opts.alias, "cea_init_thermo", char(thermolib));

    if opts.verbose
        fprintf("Initializing CEA trans: %s\n", translib);
    end
    calllib(opts.alias, "cea_init_trans", char(translib));

    ceainit = libpointer("int32Ptr", 0);
    outInit = calllib(opts.alias, "cea_is_initialized", ceainit);
    assert(outInit == "CEA_SUCCESS", "Error calling cea_is_initialized");
    assert(ceainit.Value == 1, "CEA not successfully initialized");
end