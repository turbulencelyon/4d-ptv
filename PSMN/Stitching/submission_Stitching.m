function submission_Stitching(ManipName,minframe,maxframe,dfmax,dxmax,dvmax,lmin,Session_INPUT,Session_OUTPUT)
% function submission_Stitching(ManipName,minframe,maxframe,dfmax,dxmax,dvmax,lmin,test,Session_INPUT,Session_OUTPUT)
% Function to compile to run Stitching in parallel on the PSMN
%----------------------------------------------------------------------------
%  How to use:
% Compile the file with the following command: 
% mcc -m submission_Stitching.m
% You will obtain several file, 2 are useful: run_submission_Stitching.sh and
% submission_Stitching.exe
% IMPORTANT: in run_submission_Stitching.sh, replace the line 30 by
%  eval "/My4D-PTVPath/PSMN/Tracking3D/submission_Stitching" $args
% Once you did this you can modified and then call the file
% submission_Stitching.sh in a terminal to launch the parallelisation
% --------------------------------------------------------------------------
% 2020 E. Bernard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session.input_path=Session_INPUT;
session.output_path=Session_OUTPUT;

Dfmax=str2num(dfmax)
Maxframe=str2num(maxframe)
Minframe=str2num(minframe)
Dxmax=str2num(dxmax)
Dvmax=str2num(dvmax)
Lmin=str2num(lmin)
% Test=str2num(test)
% StitchedTraj = Stitching_psmnA(session,ManipName,Minframe,Maxframe,Dfmax,Dxmax,Dvmax,Lmin,Test)
StitchedTraj = Stitching_psmnA(session,ManipName,Minframe,Maxframe,Dfmax,Dxmax,Dvmax,Lmin)
end
