function out = speed_reducer_scale(in, p, direction)
% SPEED_REDUCER_SCALE  Map design variables to/from the unit box [0,1]^7.
%   xs = speed_reducer_scale(x,  p, 'to')    x in display units
%   x  = speed_reducer_scale(xs, p, 'from')
%
% Author: AWD Labs
% Student ID: 52104479

in = in(:);
switch direction
    case 'to',   out = (in - p.lb)./(p.ub - p.lb);
    case 'from', out = p.lb + in.*(p.ub - p.lb);
    otherwise, error('speed_reducer_scale:badDirection', 'direction must be "to" or "from".');
end
end
