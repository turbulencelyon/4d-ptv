function submission_Tracking(ManipName,NbFramePerJob,minframe,maxframe,maxdist,lmin,flag_pred,npriormax,flag_conf,Session_INPUT,Session_OUTPUT)
%Used to launch the parallele calcultation in a terminal 
%----------------------------------------------------------------------------
%How to use:
%Enter your path in the session path variable (2 fields: session.input_path
% and session.output_path),
%Enter your manip name in the ManipName variable
%Compile the file with the following command: mcc -m submission_matlab.m
%You will obtain several file, 2 are usefull run_submission_matlab.sh and
%submision_matlab.exe
%IMPORTANT: in run_submission_matlab.sh, replace the line 30 by
%  eval "/Yourpath/submission_matlab" $args
%Once you did this you can modified and call the file
%submision_%tracking.sh in a terminal to launch the parallelisation



session.input_path=Session_INPUT;
session.output_path=Session_OUTPUT;


NbFramePerJob=str2num(NbFramePerJob); 
minframe=str2num(minframe);
maxframe=str2num(maxframe) ;                            
maxdist=str2num(maxdist);                          
lmin=str2num(lmin);
flag_pred=str2num(flag_pred);
npriormax=str2num(npriormax);                                                                                                                                                                          
flag_conf=str2num(flag_conf); 


[tracks,traj]=track3d_psmn(session,ManipName,NbFramePerJob,minframe,maxframe,maxdist,lmin,flag_pred,npriormax,flag_conf);

end
