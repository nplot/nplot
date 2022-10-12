function errstr = converterror(value,err,opt,prec)

% errstr = converterror(value,err,opt,prec)
%
% CONVERTERROR formats value and error into a single string
% opt:  "pm" for \pm, "textpm" for +/-, otherwise ()
% prec: precision (e.g. prec = 2 for the two first significant digits), default: 1
%
% Example: 
%   converterror(12.123456,0.004321,'pm',3)
%    12.12346 \pm 0.00432
%   converterror(12.123456,0.004321)
%    12.123(4)

% P. Steffens, 09/2012
% J. Stein, 01/2017


if nargin>2 && strcmpi(opt,'textpm')
    pmsymbol = '+/-';
    pm = true;
elseif nargin>2 && strcmpi(opt,'pm')
    pmsymbol = '\pm';
    pm = true;
else
    pm = false;
end
    
if nargin<4, prec = 1; end

if err<=0
    errstr = num2str(value);
    if pm, errstr = [errstr, ' ',pmsymbol,' 0']; end
    return;
end

% position of first significant digit in err
pos = floor(log(err)/log(10)) - prec + 1;

if round(err / 10^pos)==1, pos = pos - 1; end

pos = min([pos,0]);  % Never round before decimal point

rerr = round(err/10^pos);
rval = round(value/10^pos)  * 10^pos;
if pm, rerr = rerr * 10^pos; end

if pm
    if pos<0, fstring = ['%0.' num2str(abs(pos)) 'f']; else fstring = '%d'; end
    errstr = [num2str(rval,fstring),' ',pmsymbol,' ', num2str(rerr,fstring)];
else
    errstr = [num2str(rval) '(' num2str(rerr) ')'];
end
    



