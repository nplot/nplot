function [a4,d,names,rlist] = calcbraggangl(lattice, reflections, ki, a4minval, a4maxval)


    function ok=inlimits(angl)
        ok = angl>a4minval && angl < a4maxval;
    end

lat.a = lattice(1); lat.alpha = pi/180*lattice(4);
lat.b = lattice(2); lat.beta  = pi/180*lattice(5);
lat.c = lattice(3); lat.gamma = pi/180*lattice(6);
[~,~,B]=UBmatrix(lat,[1,0,0],[0,1,0]);

a4 = []; d= []; names = {}; rlist = {};
nolimits = nargin<5;
for n=1:size(reflections,1)
    qvec = B * reflections(n,:)';
    dval = 2*pi / sqrt(qvec'*qvec);
    if pi/dval/ki>1, continue; end
    a4val = 2*asind(pi/dval/ki);
    if nolimits || inlimits(-a4val)  % -a4
        a4 = [a4; -a4val]; d = [d; dval]; %#ok<*AGROW>
        names{end+1} = num2str(reflections(n,:));
    end
    if nolimits || inlimits(a4val)   % +a4
        a4 = [a4; a4val]; d = [d; dval];
        names{end+1} = num2str(reflections(n,:));
    end
end

% sort according to a4;
[a4,ind] = sort(a4);
d = d(ind);
names = names(ind);

% create cell array with everything
for n=1:numel(a4)
    rlist{n,1} = names{n};
    rlist{n,2} = num2str(a4(n),'% 8.2f');
    rlist{n,3} = d(n);
    rlist{n,4} = a4(n);
end

end

