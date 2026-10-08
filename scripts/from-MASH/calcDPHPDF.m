function [dt,prob] = calcDPHPDF(a,eigval,minprob)
%% [dt,prob] = calcDPHPDF(a,eigval,minprob)
%
% Calculates DPH's PMF according to the spectral decomposition in a sum of
% geometric components.
%
% a: [1-by-K] geometric components of each components
% eigval: [1-by-K] geometric constants of each component
% minprob: minimum probability mass corresponding to upper time bound in 
%          calculated PMF.
%%

% defaults
fact = 1.5;

% initializes variables
t = 1;
n = 0;
dt = [];
prob = [];
decr = false;

% calculates PMF for exponentially increasing time value
while ~isinf(t) && ~(~isempty(prob) && decr && prob(end)<minprob)
    n = n+1;
    dt = cat(1,dt,t);
    
    % calculates PMF for current time point
    y = 0;
    for k = 1:numel(eigval)
        y = y+a(k)*eigval(k)^(t-1);
    end
    prob = cat(1,prob,y);
    
    % determine whether PMF is decreasing
    if t>1 && (prob(n)-prob(n-1))<0
        decr = true;
    end
    
    % increases time point
    t = t*fact;
    
    % rounding times prevents zero disrete phase-type PMFs
    t = round(t);
    if ~isempty(dt) && t==dt(end)
        continue
    end
end
