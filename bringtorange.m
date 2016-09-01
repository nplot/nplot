function erg = bringtorange(angles)

% Problem: In a scan, the psi-value might "jump" from -180 to 180 degrees
% and thus give a discontinuous set of yang, which is a problem for the
% Voronoi cells, interpolation, etc.
% If this is the case, add 360 degrees to a part of the yang's.
erg = angles;
if max(angles(:))-min(angles(:)) > 360-4 %(4 is arbitrary, might take stdcell size **)
    % Try to find the 'hole' in the yang data
    angsort = sort(unique(angles));
    [m,i] = max(angsort(2:end)-angsort(1:end-1));
    if m>6 %(6 again arbitrary)
        % A hole exists
        hole = (angsort(i)+angsort(i+1))/2;
        % Add 360 to values left of the hole
        erg(erg<hole) = erg(erg<hole)+360;
    end
end