function [Dxinit,Dyinit,Dzinit,S,index]=DinitialComp(session,ManipName,size,track)

%% Compute the initial distance between all trajectorie of a minimum size
%Input:
%Sze= the minimum size of the trajectorie 
%track=track structure obtain after tracking3D of after stitching
%Output:
%% The output are made for the function SelectSize.m
%Dinit=the arrays of all the initial separation of all the trajectorie
%S=the structure containing all the trajectories
%index=the position of the trajectories that have the rigth size 
%------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define folder
folderout = sprintf('%sProcessed_DATA/%s/Post_Processed_Data/StructureFct/MinSize_%d/',session.output_path,ManipName,size)

if ~isfolder(folderout)
    mkdir(folderout);
end


%% Find the position of the trajectory that have the rigth size 
index=find([track.L]>size);
Ntrack=length(index);

%% Load the trajectory from a structure file that have position field (x,y,z)

Xtt=[];
Ytt=[];
Ztt=[];
fprintf('%d tracks are found for this size\n',Ntrack)
 S=struct('Id',0,'Length',0,'X',0,'Y',0,'Z',0);

for i=1:Ntrack
    elm=index(i);
    %visu_1(track,elm);
    %hold on
    X=track(elm).x;
    Y=track(elm).y;
    Z=track(elm).z;
    Ltrack=length(X);
    S(i).Id=elm;
    S(i).Length=Ltrack;
    S(i).X=X;
    S(i).Y=Y;
    S(i).Z=Z;
end 

%% Calcul des distances initiales
Dinit=[];
Dxinit=[];
Dyinit=[];
Dzinit=[];

for i=1:Ntrack
    if rem(i,1000)==0
        disp(sprintf('%d %s done',i/Ntrack*100,'%'));
    end
   Xinit1=S(i).X(1);
   Yinit1=S(i).Y(1);
   Zinit1=S(i).Z(1);
   for j=i:Ntrack
        if i~=j;
            Xinit2=S(j).X(1);
            Yinit2=S(j).Y(1);
            Zinit2=S(j).Z(1);
            Dxinit(i,j)=sqrt((Xinit1-Xinit2)^2);
            Dyinit(i,j)=sqrt((Yinit1-Yinit2)^2);
            Dzinit(i,j)=sqrt((Zinit1-Zinit2)^2);    
        end
   end
    Xinit1=[];
    Yinit1=[];
    Zinit1=[];
    Xinit2=[];
    Yinit2=[];
    Zinit2=[];
end

save(sprintf('%s/Dx%d.mat',folderout,size),'Dxinit','-v7.3')
save(sprintf('%s/Dy%d.mat',folderout,size),'Dyinit','-v7.3')
save(sprintf('%s/Dz%d.mat',folderout,size),'Dzinit','-v7.3')

save(sprintf('%s/S%d.mat',folderout,size),'S','-v7.3')
save(sprintf('%s/index%d.mat',folderout,size),'index','-v7.3')
end

