%% compare data matfiles (from different vintages) for project "kensington"

%% load toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/

%#ok<*NANVAR>
%#ok<*UNRCH>

%% clear workspace
clear variables
clear global
close all
fclose all;
clc

%% parameters
thisdir  = pwd;
thatdir  = fullfile('..', 'matdataKensington');

matfiles = dir(fullfile(thisdir, '*.mat'));

for mm = 1 : length(matfiles)
    
    fprintf('Processing %s ... ', matfiles(mm).name)
    this = load(fullfile(thisdir, matfiles(mm).name));
    that = load(fullfile(thatdir, matfiles(mm).name));

    if isequaln(this, that)
    fprintf(' OK.\n');
    else
        error('Mismatch: %s', matfiles(mm).name)

    end

end

%% finish
dockAllFigures