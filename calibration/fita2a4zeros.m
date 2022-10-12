function [newzeros,offset,calpos] = fita2a4zeros(powderdata, oldzeros, varindex, fileindex)

% PErform a fit of the a2 and a4 offsets on powder data

% P. Steffens, 07/2014

if nargin<4
    fileindex=1:numel(powderdata.fita4);
end
fileindex = fileindex(:);

x = [powderdata.scancenter(fileindex),powderdata.ki(fileindex),powderdata.dm(fileindex),powderdata.a2(fileindex)];
x = reshape(x,[],4);

%fit 
[~,opt,dpa] = funcfit(@powderpeakpos, x, reshape(powderdata.fita4(fileindex),[],1), reshape(powderdata.erra4(fileindex),[],1), [0,0], varindex);

offset.a2 = opt(1);
offset.a4 = opt(2);
newzeros.a2 = oldzeros.a2 + offset.a2 ;
newzeros.a4 = oldzeros.a4 + offset.a4 ;
if isfield(oldzeros,'a1'), newzeros.a1 = oldzeros.a1 + offset.a2/2; end
offset.da2 = dpa(1);
offset.da4 = dpa(2);

%Recalc position for all
calpos = powderpeakpos(opt,[powderdata.scancenter(:),powderdata.ki(:),powderdata.dm(:),powderdata.a2(:)]);
