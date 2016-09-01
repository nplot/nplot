function ar = polygonarea (verticesX, verticesY, faces)

% Calculate the area of polygon cells, which are given by the rows in
% "faces" - these are indices which refer to the coordinate lists verticesX/Y
% Attention: gives positive or negative results depending on order of vertices
%
% P. Steffens 05/2008

ar = zeros(size(faces,1),1);
f2 = size(faces,2);
for i=1:size(faces,1)
    ar(i)=0;
    for j = 1:f2 
        if (j==f2) || (isnan(faces(i,j+1)))
            ar(i) = ar(i) + ( verticesX(faces(i,j)) - verticesX(faces(i,1)) ) * (verticesY(faces(i,j)) + verticesY(faces(i,1)));
            break;
        else
            ar(i) = ar(i) + ( verticesX(faces(i,j)) - verticesX(faces(i,j+1)) ) * (verticesY(faces(i,j)) + verticesY(faces(i,j+1)));
        end
    end
end

ar = 1/2 * ar(:);
