function datafile(varargin)

% "datafile"-script for Matlab; calls datafile script on remote computer
%
% usage:   datafile find|show|grep filename [grep-param]
% example: datafile find 012345  (shows location of file 012345)
%          datafile find "0123*" (shows ... all matching. Note: need "")
%          datafile show 012345  (shows contents)
%          datafile grep 012345 COMND  (shows lines matching pattern COMND)

[server, username] = getoption('defaultserver','defaultuser');

system(['ssh ' username '@' server ' ./datafile ' strjoin(varargin)],'-echo');
