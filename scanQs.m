function qs = scanQs(command)

% Calculate the Qs [qh,qk,ql,en] for a scan defined by 'command'
%
% P. Steffens, 03/2008



qcent=[0,0,0,0];
dqs  =[0,0,0,0];

% Check command syntax
command = regexp(command,'(sc|bs)\s+qh\s[\s\d+-.]+dqh\s[\s\d+-.]+np\s+\d+','match');
if isempty(command), fprintf('Syntax not matched!\n'); end
command = command{1};

% 'qh'-part
qhs = regexp(command,'(?<!dqh.*)[\d.-]*','match');
for i=1:min(length(qhs),4)
    qcent(i)=str2double(qhs{i});
end

% 'dqh'-part
dqhs = regexp(command,'(?<!np.*)(?<=dqh.*)[\d.-]*','match');
for i=1:min(length(dqhs),4)
    dqs(i)=str2double(dqhs{i});
end

% 'np'-part
np = regexp(command,'(?<=np\s+)\d+','match');
np = str2double(np);

% construct scan
if strcmpi(command(1:2),'bs')
    qs = repmat(qcent,np,1) + (0:(np-1))' * dqs;
else
    qs = repmat(qcent,np,1) + (-(np-1)/2:(np-1)/2)' * dqs;
end


