function nxsfind(nxs,name)
% search the nexus file structure for a field name

% PS 06/23


function searchlevel(s,prefix)
for f = fieldnames(s)'
    f = f{1};
    if contains(f,name,'Ignorecase',true)
        fprintf('%s ',[prefix f]);
        if ischar(s.(f)), fprintf('= "%s"\n',s.(f));
        elseif isstruct(s.(f)), fprintf(': (struct)\n');
        elseif isnumeric(s.(f)), fprintf('= '); for i=1:min([3, numel(s.(f))]), fprintf('%s  ',num2str(s.(f)(i))); end, if numel(s.(f))>3, fprintf(' ...'); end, fprintf('\n');
        else, fprintf(': %s\n',class(s.(f)));
        end
    end
    if isstruct(s.(f))
        searchlevel(s.(f),[prefix f '.']);
    end
end
end

searchlevel(nxs,'.');

end
