function [fitval,optparam,dpa]=fit1d(list,dim,func,startparam,varindex,opt)

if nargin>5
    [fitval,optparam,dpa]=funcfit(func,list.coordlist(:,dim),list.valuelist(:,1),list.valuelist(:,2),startparam,varindex,opt);
else
    [fitval,optparam,dpa]=funcfit(func,list.coordlist(:,dim),list.valuelist(:,1),list.valuelist(:,2),startparam,varindex);
end