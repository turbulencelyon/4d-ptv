function submission_PairDisp(Session_INPUT,Session_OUTPUT,ManipName,SeuilInit,DeltaSeuil,DeltaInc,Epsilon,Ech,NbFrame,MinConv)
%% Compile version of PairDisp.m, usefull to launch at the PSMN 
% To compile: mcc -m submission_PairDisp.m


session.input_path=Session_INPUT;
session.output_path=Session_OUTPUT;
seuilInit=str2num(SeuilInit);
deltaSeuil=str2num(DeltaSeuil);
deltaInc=str2num(DeltaInc);
epsilon=str2num(Epsilon);
ech=str2num(Ech);
nbFrame=str2num(NbFrame);
minConv=str2num(MinConv);


Result=PairDisp(session,ManipName,seuilInit,deltaSeuil,deltaInc,epsilon,ech,minConv,nbFrame);
disp('ok');
exit;

end
