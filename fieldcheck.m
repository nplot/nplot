function erg = fieldcheck(struct, fieldnames)

% Replaces 'isfield', which does in older Matlab version not allow for cell
% arrays.

if iscell(fieldnames)
    for i=1:length(fieldnames)
        erg(i) = isfield(struct, fieldnames{i});
    end
else
    erg = isfield(struct, fieldnames);
end