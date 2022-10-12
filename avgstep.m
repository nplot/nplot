function [erg, err] = avgstep(data, step, error)

err=[];
for i=1:floor(numel(data)/step)
    if nargin > 2
        [erg(i),err(i)] = weightedmean(data(((i-1)*step+1):i*step), error(((i-1)*step+1):i*step));
    else
        erg(i) = mean(data(((i-1)*step+1):i*step));
    end
end

