function Rays_recombinaison(session,ManipName,camID)

%% Recombine the rays_camX.mat created by the function Centers2RaysParallel in one file rays.mat
%------------------------------------------------------------------------
%% Input
% session       : the directory you work in 
% ManipName     : the name of your experiment
% camID         : the list of your camera, if you have 4 cameras, it will
% be [1 2 3 4]
%
%% Output 
% A file rays.mat in the the folder Processed_DATA/MyExperiment
%-------------------------------------------------------------------------
% 2020 E. Bernard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


folderin = fullfile(session.input_path, 'Processed_DATA', ManipName);
folderout = fullfile(session.output_path, 'Processed_DATA', ManipName);

DTT=[];

fprintf("Data loading in progress")
for i=1:numel(camID)
    fileName=fullfile(session.input_path,'Processed_DATA', ManipName,['rays_cam' num2str(i) '.mat']);
    load(fileName);
    A=datacam(i).data;
    DTT(i).data=A;
    datacam=[];  
end
datacam=DTT;

fprintf("Data saving in progress...")
save(fullfile(folderout,'rays.mat'),'datacam','-v7.3')
end
