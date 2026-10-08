function [Plog, d2Plog, d1Plog] = calc_DPH_prob(T, ip, t_axis)
% [Plog, d2Plog, d1Plog] = calc_DPH_prob(T, ip, dt)
%
% Calculate the zero, first and second order derivatives of the 
% log10-probability mass function of a disrete phase-type distribution given 
% its parameters and a discrete time axis.
% Derivatives are calculated with respect to log10(x).
%
% T: [J-by-J] transition probability matrix
% ip: [1-by-J] initiation probabilities
% t_axis: [N-by-1] time axis
% Plog: [N-by-1] log10(PMF)
% d2Plog: [N-by-1] second derivative of log10(PMF) with respect to log10(t)
% d1Plog: [N-by-1] first derivative of log10(PMF) with respect to log10(t)

[a, eigval] = calcexpweight(T, ip, 0, 1);

% Define symbolic variables
syms u
x_expr = 10^(u); % Change of variable: x = 10^u, so u = log10(x)

% Symbolic construction of the log-probability
% Note: Using log() for natural log to match standard calculus derivatives
f_sym = log10(sum(a .* ((eigval').^(x_expr - 1))));

% Differentiate with respect to u (which is log10(x))
f1_sym = diff(f_sym, u, 1);
f2_sym = diff(f_sym, u, 2);

% Convert to vectorized numerical function handles
% These functions now take 'u' as an input
f_num  = matlabFunction(f_sym);
f1_num = matlabFunction(f1_sym);
f2_num = matlabFunction(f2_sym);

% To evaluate at the original points t, we must input log10(t)
u_axis = log10(t_axis);
try
    Plog = f_num(u_axis);
catch err
    Plog = [];
    d1Plog = [];
    d2Plog = [];
    return
end

try
    d1Plog = f1_num(u_axis);
catch err
    d1Plog = [];
    d2Plog = [];
    return
end

try
    d2Plog = f2_num(u_axis);
catch err
    d2Plog = [];
end