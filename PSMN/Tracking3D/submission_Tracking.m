function submission_Tracking(ManipName,NbFramePerJobMatching,minframe,maxframe,maxdist,lmin,flag_pred,npriormax,flag_conf,Session_INPUT,Session_OUTPUT)
% Function to compile to run Tracking in parallel on the PSMN
%----------------------------------------------------------------------------
%  How to use:
% Compile the file with the following command: 
% mcc -m submission_Tracking.m
% You will obtain several file, 2 are useful: run_submission_matlab.sh and
% submision_matlab.exe
% IMPORTANT: in run_submission_matlab.sh, replace the line 30 by
%  eval "/My4D-PTVPath/PSMN/Tracking3D/submission_Tracking" $args
% Once you did this you can modified and then call the file
% submision_%tracking.sh in a terminal to launch the parallelisation
% --------------------------------------------------------------------------
% 2020 E. Bernard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



session.input_path=Session_INPUT;
session.output_path=Session_OUTPUT;


NbFramePerJobMatching=str2num(NbFramePerJobMatching); 
minframe=str2num(minframe);
maxframe=str2num(maxframe) ;                            
maxdist=str2num(maxdist);                          
lmin=str2num(lmin);
flag_pred=str2num(flag_pred);
npriormax=str2num(npriormax);                                                                                                                                                                          
flag_conf=str2num(flag_conf); 


[tracks,traj]=track3d_psmn(session,ManipName,NbFramePerJobMatching,minframe,maxframe,maxdist,lmin,flag_pred,npriormax,flag_conf);

end
