function Result=PairDisp(session,ManipName,SeuilInit,DeltaSeuil,DeltaInc,Epsilon,Ech,MinConv,NbFrames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the pair dispersion for a given initial separation
% How it works:
% You will enter an initial separation, a bin size and a minimum number of couple track (for the convergence)
%and the algorithm will adapt the bin size in order to have the minimum
%number of track. The algorithm take as input the D,S and I conputed by the
%Dinitial.m function
%
%% Input:
% session        : the path where the Processed_Data is 
% ManipName      : litteraly the name of your manip
% SeuilInit      : the value of the initial separation you want to compute
% DeltaSeuil     : the bin size initial 
% DeltaInc       : the increment to reach the convergence size (same order of magnitude than SeuilInit)
% Epsilon        : the kinetic energy dissipation rate (usefull for Tstar)
% Ech            : the acquisition rate 
% Dinit          : the matrix of initial separation, output of Dinitial.m
% S              : the structure of all the tracks output of Dinitial.m
% index          : the index of all the interesseting trajectories, output of Dinitail.m
% MinConv        : the minimum number of trajectories couple you want
%
%% Output
% Result         : a structure that contain the mean of the pair
% dispersion, Tstar, DeltaSeuil (the bin size), the time and the number of
% part as a function of the time (very important to check the convergence) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Folder definition
folderout = sprintf('%sProcessed_DATA/%s/Post_Processed_Data/PairDipsersion/MinSize_%d/PairDispFinal/',session.output_path,ManipName,NbFrames)

if ~isfolder(folderout)
    mkdir(folderout);
end

D_path= sprintf('%sProcessed_DATA/%s/Post_Processed_Data/PairDipsersion/MinSize_%d/D%d.mat',session.output_path,ManipName,NbFrames,NbFrames);
S_path= sprintf('%sProcessed_DATA/%s/Post_Processed_Data/PairDipsersion/MinSize_%d/S%d.mat',session.output_path,ManipName,NbFrames,NbFrames);
I_path= sprintf('%sProcessed_DATA/%s/Post_Processed_Data/PairDipsersion/MinSize_%d/index%d.mat',session.output_path,ManipName,NbFrames,NbFrames);

load(D_path);
load(S_path);
load(I_path);
%% Convergence limit
A=find(Dinit<SeuilInit+DeltaSeuil & Dinit>SeuilInit);
L=length(A);

while L<MinConv
   DeltaSeuil=DeltaSeuil+DeltaInc
   
   A=find(Dinit<SeuilInit+DeltaSeuil & Dinit>SeuilInit);
   L=length(A)
end

%% Call the function that will compute the distance between two tracks
tic;
Delta=SelectSize(S,Dinit,index,SeuilInit,SeuilInit+DeltaSeuil);
toc;
Delta(Delta==0)=nan;
Dm=nanmean(Delta);

for k=1:length(Dm)
    B=Delta(:,k);
    NbPart(k)=numel(find(~isnan(B)==1));
end


%% Time definition 
Time=linspace(0,length(Dm)/Ech,length(Dm));
Tstar=((((2*(SeuilInit)+DeltaSeuil)/2)*10^-3)^2/(Epsilon))^(1/3);


%% Result saving
Result.Dinit=SeuilInit;
Result.BinSize=DeltaSeuil;
Result.Tstar=Tstar;
Result.Time=Time;
Result.MeanD=Dm;
Result.NbPart=NbPart;

save(sprintf('%s/pair_disp_%d.mat',folderout,SeuilInit),'Result','-v7.3');
end