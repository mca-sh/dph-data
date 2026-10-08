function [dt,prob] = calcCPHPDF(a,eigval,minprob)
%% [dt,prob] = calcCPHPDF(a,eigval,minprob)
%
% Calculates CPH's PDF according to the spectral decomposition in a sum of
% exponential components.
%
% a: [1-by-K] exponential coefficients of each component
% eigval: [1-by-K] exponential constants of each component
% minprob: minimum probability corresponding to upper time bound in 
%          calculated PDF.
%%

% defaults
fact = 1.5;

% initializes variables
t = 0;
n = 0;
dt = [];
prob = [];
decr = false;

% calculates PDF for exponentially increasing time value
while ~isinf(t) && ~(~isempty(prob) && decr && prob(end)<minprob)
    n = n+1;
    dt = cat(1,dt,t);
    
    % calculates PDF for current time point
    y = 0;
    for k = 1:numel(eigval)
        y = y+a(k)*exp(eigval(k)*t);
    end
    
    % dtermines whether PDF is decreasing
    prob = cat(1,prob,y);
    if t>0 && (prob(n)-prob(n-1))<0
        decr = true;
    end
    
    % increases time point
    if t==0
        t = 1;
    else
        t = t*fact;
    end
end
