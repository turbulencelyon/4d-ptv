function submission_CenterFinding2D(ManipName,CamNum,FirstFrame,Nframes,Th,Sz,Session_INPUT,Session_OUTPUT)

%%% Define the argument use by the function CenterFunction2D to run it
%%% in a terminal at PSMN 
%--------------------------------------------------------------------------------
%%%     ManipName       : Name of the folder experiment
%%%     CamNum          : number of the camera studied
%%%     FirstFrame      : first frame studied
%%%     Th              : threshold
%%%     Sz              : typical size of the particles
%%%     Session_INPUT   : The path of the DATA directory  
%%%     Session_OUTPUT  : The path of the PROCESSED DATA directory 

%--------------------------------------------------------------------------------
%   How to run: 
%   Compile in matlab using: mcc -m submission_CenterFinding2D.m
%   Modify the path to submission_CenterFinding2D
%   Launch it in a terminal using: sh run_submission_CenterFinding2D.sh $MCRROOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


session.input_path=Session_INPUT;
session.output_path=Session_OUTPUT;
camNum=str2num(CamNum);
firstframe=str2num(FirstFrame);
nframes=str2num(Nframes);
th=str2num(Th);
sz=str2num(Sz);
Test=0;



fprintf('Start CenterFinding2D computation... \n');
CC = CenterFinding2D(session,ManipName,camNum,firstframe,nframes,th,sz,Test);
end
