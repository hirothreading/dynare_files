function tests = test_mod_files
tests = functiontests(localfunctions);
end

function testEstimationAndSimulationShareModel(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fileparts(testFolder);
estimateText = fileread(fullfile(rootFolder, 'Estimation', 'estimate.mod'));
simulateText = fileread(fullfile(rootFolder, 'Simulation', 'simulate.mod'));
includeLine = '@#include "../Model/bll_model.inc"';
verifySubstring(testCase, estimateText, includeLine);
verifySubstring(testCase, simulateText, includeLine);
end

function testEstimatedParametersAreNotCalibratedInSteadyState(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fileparts(testFolder);
steadyStateText = fileread(fullfile(rootFolder, 'Estimation', 'estimate_steadystate.m'));
verifyFalse(testCase, contains(steadyStateText, 'cal_w ='));
verifyFalse(testCase, contains(steadyStateText, 'zet ='));
verifyFalse(testCase, contains(steadyStateText, 'sig_w ='));
end

function testWageShockUsesPublishedReplicationUnits(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fileparts(testFolder);
modelText = fileread(fullfile(rootFolder, 'Model', 'bll_model.inc'));
verifySubstring(testCase, modelText, '+ m_w;');
verifyFalse(testCase, contains(modelText, 'kap_w'));
end

function testPaperProfileContainsAppendixMarkupNormalization(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fileparts(testFolder);
modelText = fileread(fullfile(rootFolder, 'Model', 'bll_model.inc'));
verifySubstring(testCase, modelText, '@#if PAPER_SPECIFICATION == 1');
verifySubstring(testCase, modelText, '))*m_p;');
verifySubstring(testCase, modelText, '))*m_w;');
end

function testGovernmentShareIsNotEstimated(testCase)
testFolder = fileparts(mfilename('fullpath'));
rootFolder = fileparts(testFolder);
estimateText = fileread(fullfile(rootFolder, 'Estimation', 'estimate.mod'));
estimatedBlock = extractBetween(estimateText, 'estimated_params;', 'end;');
verifyFalse(testCase, contains(estimatedBlock{1}, 'psi,'));
end
