function [eqset, eqvalues, deviate, notall] = checkconstants(datalist)

% Check if constant values of the lists are consistent
% datalist is cell array
% eqset, eqvalues : names and values of constants that appear in ALL datalist{i}
% deviate         : names of constants that differ too much
% notall          : names of constants that do not appear in all
%
% P. Steffens, 08/2008



eqset = {};
eqvalues = {};
deviate = {};
notall = {};
constlist = {};

maxdeviate = getoption('maxdeviate');

% Get the names of all constants that appear
for i = 1:length(datalist)
    if isfield(datalist{i},'constants'), constlist = {constlist{:}, datalist{i}.constants{:}}; end
end
constlist = unique(constlist);

for j=1:length(constlist) % for every constant
    all = true;
    values = [];
    for i = 1:length(datalist)
        if isfield(datalist{i},constlist{j})
            values = [values; datalist{i}.(constlist{j})(:)']; %collect values 
        else
            all = false; % This constant does not appear in all data sets
            notall = {notall{:}, constlist{j}};
        end
    end
    if  ~isempty(values) && any(max(values,[],1) - min(values,[],1) > maxdeviate.(constlist{j}))
        % Deviation is too large
        deviate = {deviate{:}, constlist{j}};
    elseif all
        % if deviation ok AND existent in all data sets:
        eqset = {eqset{:}, constlist{j}};
        eqvalues = {eqvalues{:}, mean(values,1)};
    end            
end

