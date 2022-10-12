function showfunction(fname)

% Show function information

if isempty(fname)
    fprintf('Usage: showfunction [name], where [name] is a valid function .m-file.\n');
    return;
end
    
try
    eval([' [val, paramnames, paramnum, description] = ' fname '([],[],''info'');']);
    fprintf('Name:         %s\nDescription:  %s\n', fname, description);
    fprintf('Parameter names and order:\n');
    for i=1:paramnum, fprintf('   %2d: %s\n',i,paramnames{i}); end
    fprintf('Filename: %s\n',which(fname))
catch
    fprintf('Error: function %s is not found or does not have the right syntax.\n',fname);
end
    
