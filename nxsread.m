function data = nxsread(filename)

    function dat = readgroup(gindex)

        dat = [];
        group = fileinfo;
        ind = 0;
        while ind<numel(gindex)
            ind = ind+1;
            group = group.Groups(gindex(ind));
        end
        gname = group.Name;

        gid = H5G.open(fid,gname); 
        
        % cycle through Datasets
        for i=1:length(group.Datasets)
            did = H5D.open(gid,group.Datasets(i).Name);
            val = H5D.read(did);
            val = squeeze(val); % remove extra dimensions
            if ischar(val) && size(val,1)>1 && size(val,2)==1, val=val'; end % transpose to make simple string
            dat.(group.Datasets(i).Name) = val;
            H5D.close(did);
        end
        
        % cycle through subgroups
        for i=1:length(group.Groups)
            thisgname = strsplit(group.Groups(i).Name,'/');  thisgname = thisgname{end};
            dat.(thisgname) = readgroup([gindex,i]);
        end

        H5G.close(gid);
    end

    if ~exist(filename,'file'), warning('File does not exist.'); data=[]; return; end

    fid = H5F.open(filename);

    fileinfo = h5info(filename);

    data = readgroup(1);

    H5F.close(fid);

    for ii=1:length(fileinfo.Attributes)
        data.(fileinfo.Attributes(ii).Name) = fileinfo.Attributes(ii).Value;
    end

    
end


