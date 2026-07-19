%% Prepare a positive-definite Metropolis proposal from the computed mode
% The posterior mode for gamma_pi lies on the determinacy boundary. Dynare's
% numerical Hessian can therefore have one negative eigenvalue even when the
% remaining curvature is well behaved. This script changes only the proposal
% Hessian; it does not change the mode, likelihood, priors, or posterior.

clearvars;
sourceFile = fullfile('estimate', 'Output', 'estimate_mode.mat');
destinationFolder = 'mode';
destinationFile = fullfile(destinationFolder, 'estimate_mode.mat');

if ~isfile(sourceFile)
    error('BLL:MissingMode', ...
        'Run estimate.mod with MH_REPLIC=0 before preparing the mode file.');
end
modeResult = load(sourceFile, 'xparam1', 'hh', 'parameter_names', 'fval');
if any(~isfinite(modeResult.xparam1)) || any(~isfinite(modeResult.hh), 'all')
    error('BLL:InvalidMode', ...
        'The computed mode or Hessian contains non-finite values.');
end

symmetrizedHessian = (modeResult.hh+modeResult.hh')/2;
[eigenvectors, eigenvalues] = eig(symmetrizedHessian, 'vector');
positiveEigenvalues = eigenvalues(eigenvalues > 0);
if isempty(positiveEigenvalues)
    error('BLL:InvalidHessian', 'The mode Hessian has no positive eigenvalues.');
end

sortedCurvature = sort(positiveEigenvalues);
replacementCurvature = sortedCurvature(max(1, ceil(0.1*numel(sortedCurvature))));
regularizedEigenvalues = max(eigenvalues, replacementCurvature);
hh = eigenvectors*diag(regularizedEigenvalues)*eigenvectors';
hh = (hh+hh')/2;
xparam1 = modeResult.xparam1;
parameter_names = modeResult.parameter_names;
fval = modeResult.fval;

if ~isfolder(destinationFolder)
    mkdir(destinationFolder);
end
save(destinationFile, 'xparam1', 'hh', 'parameter_names', 'fval');
fprintf(['Prepared %s\nReplaced %d non-positive eigenvalue(s); ' ...
    'minimum proposal curvature is %.6g.\n'], destinationFile, ...
    sum(eigenvalues <= 0), min(regularizedEigenvalues));
