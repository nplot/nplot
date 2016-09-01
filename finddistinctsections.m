function [nsecs,section] = finddistinctsections(faces)

% function [section,nsecs] = finddistinctsections(faces)
%
% Find distinct (unconnected) sections. 
% nsecs :   number of distinct sections
% section : for each face (=line in faces), a number (indicating number of section)

% P. Steffens 01/2012

nfac = size(faces,1); fcol = size(faces,2);
section = zeros(nfac,1);
nextsec = 0;
for f = 1 : nfac
    if section(f)==0
        nextsec = nextsec + 1;
        section(f) = nextsec;
    end
    for v = faces(f, isfinite(faces(f,:)))
        fi = find(faces==v);
        for fv = mod(fi'-1,nfac)+1 % floor((fi'-1)/fcol)+1; % faces containing this vertex v
            if section(fv)>0 && section(fv)~=section(f)
                section(section==section(fv))=section(f);
            else
                section(fv) = section(f);
            end
        end
        % all faces containing vertex v are now numbered like face f
        faces(fi) = nan; % set to nan to avoid further iterations on it
    end
end

[un,u,section] = unique(section); % renumber
nsecs = numel(un);  % number of distinct sections        