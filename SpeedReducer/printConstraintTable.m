function printConstraintTable(x, params)
% printConstraintTable  Print every constraint: value, limit, utilisation, status.
% Author: AWD Labs
%
% Status is ACTIVE if the utilisation is within params.activeTol of 100 %
% (for example activeTol = 1e-3 means 99.9 % or more counts as active).
% A constraint is VIOLATED if it is over its limit by more than params.feasTol.

r = speedReducerAnalysis(x, params);
c = constraintFun(x, params);
utilisation = (c + 1) * 100;      % percent

fprintf('%-28s %10s %10s %-5s %12s   %s\n', ...
        'Constraint', 'Value', 'Limit', 'Unit', 'Utilisation', 'Status');
for k = 1:numel(c)
    if c(k) > params.feasTol
        status = 'VIOLATED';
    elseif c(k) >= -params.activeTol
        status = 'ACTIVE';
    else
        status = 'INACTIVE';
    end
    fprintf('%-28s %10.4f %10.4f %-5s %10.2f %%   %s\n', ...
            params.conNames{k}, ...
            r.conValue(k) * params.conToDisp(k), ...
            r.conLimit(k) * params.conToDisp(k), ...
            params.conUnits{k}, utilisation(k), status);
end
end
