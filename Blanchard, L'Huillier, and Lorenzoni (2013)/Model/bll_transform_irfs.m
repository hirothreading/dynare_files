function structuralIRF = bll_transform_irfs(reducedIRF, rho, sigU, sigNu)
%BLL_TRANSFORM_IRFS Map EFI innovations into structural technology shocks.
%   reducedIRF is T-by-2 and contains responses to e_1 and e_2. The output
%   is T-by-3 for permanent technology, transitory technology, and noise.

arguments
    reducedIRF (:,2) double {mustBeFinite}
    rho (1,1) double
    sigU (1,1) double
    sigNu (1,1) double
end

[~, ~, C, filter] = bll_signal_extraction(rho, sigU, sigNu);
A = filter.A;
F = filter.F;
K = filter.K;
T = size(reducedIRF, 1);
structuralIRF = zeros(T, 3);

for shock = 1:3
    if shock == 1
        stateShock = [(1-rho)*sigU; 0; 0];
        signalShock = [0; 0];
    elseif shock == 2
        stateShock = [0; 0; sqrt(rho)*sigU];
        signalShock = [0; 0];
    else
        stateShock = [0; 0; 0];
        signalShock = [0; sigNu];
    end

    state = zeros(3, T);
    estimatedState = zeros(3, T);
    signal = zeros(2, T);
    innovation = zeros(2, T);
    orthogonalInnovation = zeros(2, T);

    state(:,1) = stateShock;
    signal(:,1) = F*state(:,1) + signalShock;
    estimatedState(:,1) = K*signal(:,1);
    innovation(:,1) = signal(:,1);
    orthogonalInnovation(:,1) = C\innovation(:,1);
    structuralIRF(1,shock) = reducedIRF(1,:)*orthogonalInnovation(:,1);

    for period = 2:T
        state(:,period) = A*state(:,period-1);
        signal(:,period) = F*state(:,period);
        forecast = F*A*estimatedState(:,period-1);
        innovation(:,period) = signal(:,period)-forecast;
        estimatedState(:,period) = A*estimatedState(:,period-1) ...
            + K*innovation(:,period);
        orthogonalInnovation(:,period) = C\innovation(:,period);
        structuralIRF(period,shock) = sum( ...
            reducedIRF(1:period,:).*flipud(orthogonalInnovation(:,1:period)'), ...
            'all');
    end
end
end
