function [ismin,ismax] = islocalextr(dat)
% [ismin,ismax] = islocalextr(dat)
%
% Identifies local maxima in input series.
%
% dat: data vector
% ismin: logical vector containing trues at local minima
% ismax: logical vector containing trues at local maxima

% Size must be initialized to account for empty dat vector
ismax = false(size(dat));
ismin = false(size(dat));

% This works also for empty dat vector
ismax(2:end-1) = dat(2:end-1)>dat(3:end) & dat(2:end-1)>dat(1:end-2);
ismin(2:end-1) = dat(2:end-1)<dat(3:end) & dat(2:end-1)<dat(1:end-2);