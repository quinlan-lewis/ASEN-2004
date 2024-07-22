
function [value,isterminal,direction] = stopping_point(t,X)
% Locate the time when height passes through zero in a decreasing direction
%and stop integration.
value = X(2); % detect when z = X(3) < 0
isterminal = 1; % stop the integration
direction = -1; % negative direction
end

