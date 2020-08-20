function Submission_Matlab_Stitching(ManipName,minframe,maxframe,dfmax,dxmax,dvmax,lmin,Session_INPUT,Session_OUTPUT)


session.input_path=Session_INPUT;
session.output_path=Session_OUTPUT;

Dfmax=str2num(dfmax)
Maxframe=str2num(maxframe)
Minframe=str2num(minframe)
Dxmax=str2num(dxmax)
Dvmax=str2num(dvmax)
Lmin=str2num(lmin)
StitchedTraj = Stitching_psmnA(session,ManipName,Minframe,Maxframe,Dfmax,Dxmax,Dvmax,Lmin)
end
