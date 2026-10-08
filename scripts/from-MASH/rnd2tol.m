function val = rnd2tol(val0,tol)
%% val = rnd2tol(val0,tol)
%
% Rounds up value to closest integer tolerating a certain minimal deviation
% 
% val0: value to round up
% tol: tolerated deviation
% val: rounded value
%%
if tol==0
    val = val0;
    return
end
val = tol*round(real(val0)/tol)+1i*tol*round(imag(val0)/tol);