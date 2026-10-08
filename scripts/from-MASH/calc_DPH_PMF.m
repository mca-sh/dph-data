function prob = calc_DPH_PMF(T,t,ip,dt)
ndt = size(dt,2);
prob = zeros(1,ndt);
for n = 1:ndt
    prob(n) = ip*(T^(dt(n)-1))*t;
end