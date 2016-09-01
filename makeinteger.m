function newvec = makeinteger(vec,prec)

% try something to make vector components integer...

vec = round(vec/prec)*prec;
if all(vec==0), newvec=vec; return; end
nz = vec~=0;
scales = abs(1./vec(nz));
ps = scales(1):scales(1):10/max(abs(vec));
for d=2:sum(nz)
    ps2 = scales(d):scales(d):10/max(abs(vec));
    in = intersect(ps,ps2);
    if isempty(in), ps=max([min(ps),min(ps2)]); else ps=in; end
end
ps = min(ps);
if isempty(ps),ps=1; end
newvec = round(ps*vec/prec)*prec;


       