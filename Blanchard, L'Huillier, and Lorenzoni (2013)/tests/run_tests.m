%% Run the fast MATLAB regression suite
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fileparts(testFolder);
addpath(fullfile(rootFolder, 'Model'));
addpath(fullfile(rootFolder, 'Estimation'));
results = runtests(testFolder, 'IncludeSubfolders', false);
assertSuccess(results);
disp(table(results));
