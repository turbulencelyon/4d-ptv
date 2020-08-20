function Rays_recombinaison(session,ManipName,camID)

%% Recombine the rays_camX.mat created by the function Centers2RaysParallel in one file rays.mat
%
%% Input
% session       : the directory you work in 
% ManipName     : the name of your experiment
% camID         : the list of your camera, if you have 4 cameras, it will
% be [1 2 3 4]
%
%% Output 
% A file rays.mat in the the folder Processed_DATA/MyExperiment
%%


folderin = sprintf("%sProcessed_DATA/%s/",session.input_path,ManipName);
folderout = sprintf("%sProcessed_DATA/%s/",session.output_path,ManipName);

DTT=[];

for i=1:numel(camID)
    fileName=sprintf('%sProcessed_DATA/%s/rays_cam%d.mat',session.input_path,ManipName,i);
    load(fileName);
    A=datacam(i).data;
    DTT(i).data=A;
    datacam=[];  
end
datacam=DTT;

save(sprintf('%s/rays.mat',folderout),'datacam','-v7.3')
end
