function Kinf = HinfRstability_regulation(A,B1,B2,C1,D11,D12, alpha, r, theta)


% Sizes
nx = size(A,1);
nu = size(B1,2);
nd = size(B2,2);
ny = size(C1,1);

% YALMIP decision variables
X = sdpvar(nx, nx, 'symmetric'); % Symmetric positive definite matrix X
Y = sdpvar(nu, nx);              % Matrix Y
gamma = sdpvar(1,1);

% Constraints
LMI1 = X >= 0;

LMI2 = [(A*X+B1*Y)+(A*X+B1*Y)' B2 (C1*X+D12*Y)';
        B2' -gamma*eye(nd) D11';
        (C1*X+D12*Y) D11 -gamma*eye(ny)] <= 0;

LMI3 = 2 * alpha * X + (A * X + B1 * Y) + (A * X + B1 * Y)' <= 0;

LMI4 = [ -r * X,           (A * X + B1 * Y);
         (A * X + B1 * Y)', -r * X       ] <= 0;

sinTheta = sin(theta);
cosTheta = cos(theta);
LMI5 = [ ((A * X + B1 * Y) + (A * X + B1 * Y)') * sinTheta, ((A * X + B1 * Y) - (A * X + B1 * Y)') * cosTheta;
         ((A * X + B1 * Y)' - (A * X + B1 * Y)) * cosTheta, ((A * X + B1 * Y)' + (A * X + B1 * Y)) * sinTheta ] <= 0;

con = [LMI1 LMI2 LMI3 LMI4 LMI5];

% Objective
obj = gamma;

% Options
opt = sdpsettings('solver', 'lmilab', 'verbose', 0);

sol = optimize(con, obj, opt);

% Check solver status
if sol.problem == 0
    X_sol = value(X);
    Y_sol = value(Y);
    
    % Compute controller gain K
    Kinf = Y_sol / X_sol;
else
    error('Failed to solve the LMIs. Solver message: %s', sol.info);
end

end

