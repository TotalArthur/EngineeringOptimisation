function f = objectiveFun(x, params)
% objectiveFun  Objective for fmincon: total material volume of the reducer.
% Author: AWD Labs
%
% x is the design vector in brief units [b m z l1 l2 d1 d2] (cm, mm, -, cm...).
% The volume is returned in cm^3 divided by params.objScale. Dividing by
% 1000 just keeps the objective near 1 to 5, which helps the solver. Multiply
% by params.objScale to get cm^3 back.

r = speedReducerAnalysis(x, params);
f = r.V_total * params.toCm3 / params.objScale;
end
