function [BM, FA, CM, diagnostics] = bll_signal_extraction(rho, sigU, sigNu)
%BLL_SIGNAL_EXTRACTION Solve the steady-state consumer Kalman filter.
%   The returned matrices are the scalar-compatible replacement for the
%   matrix expressions in the authors' pre-Dynare 4.5 model block.

arguments
    rho (1,1) double {mustBeFinite, mustBeGreaterThan(rho,0), mustBeLessThan(rho,1)}
    sigU (1,1) double {mustBeFinite, mustBePositive}
    sigNu (1,1) double {mustBeFinite, mustBePositive}
end

A = [1+rho, -rho, 0; 1, 0, 0; 0, 0, rho];
F = [1, 0, 1; 1, 0, 0];
L1 = diag([(1-rho)*sigU, 0, sqrt(rho)*sigU]);
L2 = diag([0, sigNu]);
S1 = L1*L1';
S2 = L2*L2';

P = 0.01*eye(3);
tolerance = 1e-13;
maxIterations = 100000;
converged = false;

innovationCovariance = F*P*F' + S2;
K = (P*F')/innovationCovariance;
for iteration = 1:maxIterations
    Pnext = A*(P-K*F*P)*A' + S1;
    Pnext = (Pnext+Pnext')/2;
    innovationCovariance = F*Pnext*F' + S2;
    Knext = (Pnext*F')/innovationCovariance;
    convergenceResidual = norm(Knext-K, Inf);
    if convergenceResidual < tolerance
        P = Pnext;
        K = Knext;
        converged = true;
        break
    end
    P = Pnext;
    K = Knext;
end

if ~converged
    error('BLL:KalmanNoConvergence', ...
        'The signal-extraction Riccati iteration did not converge.');
end
if any(~isfinite(P), 'all') || any(~isfinite(K), 'all')
    error('BLL:InvalidKalmanSolution', ...
        'The signal-extraction solution contains non-finite values.');
end

innovationCovariance = (F*P*F'+S2);
innovationCovariance = (innovationCovariance+innovationCovariance')/2;
[eigenvectors, eigenvalues] = eig(innovationCovariance, 'vector');
scale = max(1, max(abs(eigenvalues)));
if min(eigenvalues) < -1e-12*scale
    error('BLL:InvalidInnovationCovariance', ...
        'The innovation covariance matrix is not positive semidefinite.');
end
eigenvalues = max(eigenvalues, 0);
CM = eigenvectors*diag(sqrt(eigenvalues))*eigenvectors';
CM = (CM+CM')/2;

BM = K*CM;
FA = F*A;

diagnostics = struct( ...
    'A', A, 'F', F, 'K', K, 'P', P, ...
    'innovationCovariance', innovationCovariance, ...
    'iterations', iteration, 'residual', convergenceResidual);
end
