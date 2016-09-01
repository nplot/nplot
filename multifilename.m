function filelist = multifilename(file)

% Returns cell array of full path filenames

filelist = [];
multifile=[];

%----- Find if there's multiple scans

if ~isempty(findstr(file,'['))
    
    checkstring = file( (findstr(file,'[')) : (findstr(file,']')) );

    % let Matlab evaluate string
    try
        specs = eval(checkstring);
    catch
        fprintf('String cannot be evaluated. Check input!\n'); return;
    end
    
    % check for leading zeroes in all elements
    z1 = regexp(checkstring,'\d+'); % starting ind. of all elements
    z2 = regexp(checkstring,'(?=0*)[1-9]\d*');  % -"- without leading zeroes
    zeroadd = min(z2-z1);

    filefront = file(1:findstr(file,'[')-1);

    fileback = file(findstr(file,']')+1:end);

    for i = 1:length(specs)

        multifile = [multifile; filefront num2str(specs(i),['%0' num2str(ceil(log10(max(specs)+1)))+zeroadd 'd']) fileback];  %ensure same number of digits in string

    end

else

    multifile = file;

end



%----- Put the file numbers into an array of cells

filelist = cell(size(multifile(:,1)));

nf = length(filelist);



for j = 1:nf

    filelist{j} = multifile(j,:);

end

