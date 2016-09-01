function [newvertices, newfaces, vassign] = optimizepatch(vertices, faces)

% Makes vertex list unique and removes unused vertices 
% (correct assignment in faces accordingly)
% vassign is assignment between old and new vertexlist (newvertices = vertices(vassign,:))

% P. Steffens, 05/2008 - 08/2014


faces(faces==0) = nan;
newfaces = faces;

%reduce vertex list to unique values (within precision 1E-10)
[newvertices,m1,index]=unique(round(vertices*1E10)/1E10,'rows'); 
%...and replace indices in face list correspondingly
newfaces(isfinite(faces))=index(faces(isfinite(faces))); 

% Check for doubles ??? ** not yet... (do outside)

% any empty columns?
newfaces = newfaces(:, any(~isnan(newfaces),1));

% Try to delete unused vertices...
usedvert = newfaces(isfinite(newfaces));
[usedvert,~,index2]=unique(usedvert); 
newvertices = newvertices(usedvert,:);
newfaces(isfinite(newfaces)) = index2;


vassign = m1(usedvert);
