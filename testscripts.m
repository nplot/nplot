% d=nxsread('036806.nxs'); % IN20
d=nxsread('044853.nxs'); % THALES
% d=nxsread('047918.nxs'); % IN8
% d=nxsread('072803.nxs'); % IN12
% d=nxsread('061580.nxs'); % IN22


nxsfind(d,'A4')
nxsgetvar(d,'A4')
t=nxsTAS(d)
nxsprintscan(d)