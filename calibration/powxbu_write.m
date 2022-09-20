function powxbu_write(filename, ki, d, a4, da4, np, ti, comment, direction)

% write xbu file for a4 powder scans, fot given reflection list

% P. Steffens, 07/2014

a4 = a4(:);
[~, a4, ki, d, da4, np, ti] = makesamesize(a4, ki, d, da4, np, ti);
if nargin<8, comment = []; end;
lastki = nan;

try
    file = fopen(filename,'w');
catch
    fprintf('Error: could not open file %s for writing.\n',filename);
    return;
end

% sort list according to ki and a4
[~,order] = sortrows([ki,a4]);
if nargin>8 && direction<0, order = order(end:-1:1); end
a4 = a4(order); ki = ki(order); d = d(order); da4 = da4(order); np = np(order); ti = ti(order);
% detect sign changes in a4
changeind = sign(a4(2:end))~=sign(a4(1:end-1));
changeind = find(changeind); % indices of those a4 ofter which sign changes



fprintf(file,'! %s\n! Created on %sh%s by powxbu\n\n',comment,datestr(now,'yyyy-mm-dd at HH'),datestr(now,'MM'));
fprintf(file,'se title Powder\n\n');
fprintf(file,'se fx 1\n\n');
fprintf(file,'! Drive to starting point\n');
fprintf(file,'dr a6 5\n');
fprintf(file,'dr a4 %7.2f\n', a4(1) - (np(1)-1)/2*da4(1));
fprintf(file,'dr a6 0\n');

for n = 1:numel(a4)
    if ki(n)~=lastki % drive ki
        fprintf(file,'\n! drive KI\ndr ki %7.3f\n\n',ki(n));
        lastki = ki(n);
    end
    % a4-scan
    fprintf(file,'se as %7.4f\n',d(n));
    fprintf(file,'sc a4 %7.2f da4 %5.2f np %3d ti %3d\n',a4(n), da4(n), np(n), ti(n));
    if ismember(n,changeind) % a4 goes to other side
        if a4(n+1)>0, fprintf(file,'\n! go to positive a4 (turn a6)\n');
        else fprintf(file,'\n! go to negative a4 (turn a6)\n'); end
        fprintf(file,'dr a6 5\n');
        fprintf(file,'dr a4 %7.2f\n', a4(n+1) - (np(n+1)-1)/2*da4(n+1));
        fprintf(file,'dr a6 0\n\n');
    end
end

fprintf(file,'\nse fx 2\n\n! Do not forget to reset TITLE and lattice constant AS before continuing.\n');

fclose(file);
        
