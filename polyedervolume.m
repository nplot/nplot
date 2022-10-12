function vol = polyedervolume( vertices, cells)

% Volume of N-dim polyeders defined by cells
%
% P. Steffens, 08/2008

ncells = size(cells,1);
ndims = size(vertices,2);

if ndims == 2
    try
    vol = abs ( polygonarea (vertices(:,1), vertices(:,2), cells) );
    catch
        1
    end
    % in 2D, use simpler algorithm
else
    vol = zeros(ncells,1);
    for i = 1:ncells
        [K,vol(i)] = convhulln( vertices( cells(i, isfinite(cells(i,:))), :) );
    end
end
