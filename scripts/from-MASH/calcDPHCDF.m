function cdf = calcDPHCDF(T,ip,dt,tol)
% cdf = calcDPHCDF(T,ip,dt,tol)
%
% Calculate and returns DPH's CDF.
%
% T: [D-by-D] transition probability matrix between degenerate states
% ip: [1-by-D] starting probbailities
% dt: [1-by-ndt] time axis
% tol: precision on PDF values
% cdf: [1-by-ndt] CDF values

ndt = size(dt,2);
cdf = zeros(1,ndt);
D = size(T,1);
for n = 1:ndt
    cdf(n) = 1-rnd2tol(ip*(T^dt(n))*ones(D,1),tol);
end