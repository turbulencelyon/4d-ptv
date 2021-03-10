function submission_Centers2Rays(kcam,Calib,ManipName,FirstFrame,Session_INPUT,Session_OUTPUT)
%   Define here the argument used by the function Center2Rays.m
%----------------------------------------------------------------------------------------
% Parameters:
%   session     : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
%   ManipName   : Name of the folder experiment
%   Calib       : Path of the calibration file
%   camID       : List of camera numbers. ex: [1,2,3] if you have 3 cameras numbered 1,2,3 respectively.
% ------------------------------------------------------------------------------------------
% 2020 E. Bernard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
session.input_path=Session_INPUT;
session.output_path=Session_OUTPUT;

kcam=str2num(kcam);
firstframe=str2num(FirstFrame)
[P,V]=Centers2RaysParallel(session,ManipName,Calib,kcam,firstframe);
toc
end 