function [A, b] = constraintmatrix (constring)

% Interpretes a string containing linear equation in the form " p1+p2=3; p3=5; etc."
% and return a Matrix A and a vector b to write the constraints in the form A * p = b

% P. Steffens, 10/2008


A = []; b = [];
constring = upper(strtrim(constring));
paramnum = 0;

[st,en] = regexp(constring, 'P\d+'); %find all Px entries in calcstring
for i=numel(st):-1:1
    % Convert constring in matlab readable string (PALx --> PAL(:,x))
    pnum = constring(st(i)+1:en(i));
    constring = [ constring(1:st(i)) '(:,' pnum ')' constring(en(i)+1:end) ];
    paramnum = max(paramnum, str2double(pnum));
end
if ~strcmp(constring(end),';'), constring(end+1) = ';'; end   % Add ; at the end

% Break into lines
st = findstr(constring, ';'); 
nlines = numel(st);
for i=1:nlines
    zeile{i} = constring(1:st(i));
    constring = constring(st(i)+1:end);
    st = st - st(i);
end

% Loop over lines
for nl = 1:nlines
    thisline = zeile{nl};
    % Bring everything on left side of equation
    st = findstr(thisline,'=');
    try thisline = [thisline(1:st-1) '-(' thisline(st+1:end-1) ')'];
    catch, fprintf('Error: wrong format in constraint, line %d\n', nl); return; end
    
    % Get constants
    P = zeros(1,paramnum); %#ok<NASGU>
    try  eval(['b(nl) = -(' thisline ') ; ']); 
    catch, fprintf('Error: could not interpret constraint, line %d\n', nl); A = []; b = []; return; end
    
    % Get coefficients
    P = eye(paramnum); %#ok<NASGU>
    try eval(['A(nl,:) = (' thisline ')'' + b(nl) ; ']);
    catch, fprintf('Error: could not interpret constraint, line %d\n', nl); A = []; b = []; return; end
    
end

b=b';
