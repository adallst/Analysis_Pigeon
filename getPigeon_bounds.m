function vals = getPigeon_bounds(steps, ndt)
% function vals = getPigeon_bounds(steps, ndt)
%
% Returns:
%   vals(1:3) ... [bound boundLo boundHi]
%   vals(4)   ... decision time, does not include ndt
if isempty(steps)
    vals = [0 0 0 0];
elseif length(steps) == 1
    vals = [steps(1)./2 0 steps(1) 1];
else
    ndt = min(ndt, length(steps)-2);
    steps = steps(1:end-ndt);
    vals = [mean(steps(end-1:end)), steps(end-1), steps(end) length(steps)];
end