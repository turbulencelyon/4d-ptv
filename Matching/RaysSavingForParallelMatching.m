function RaysSavingForParallelMatching(session,ManipName,camID,nbFramePerJobMatching)
%%% Save data from rays.mat into multiple small .dat to run later parallele matching
%----------------------------------------------------------------------------
%%% Parameters : 
%%%     session               : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
%%%     ManipName             : Name of the folder experiment
%%%     camID                 : List of cameras number
%%%     NbFramePerJobMatching : number of frames per job. Pay attention, has to be
%%%     chosen as a function of processing time of one picture, in order to
%%%     that each job runs for 10 min (PSMN requirements).
%------------------------------------------------------------------------------
% 2020-2021 D. Dumont
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Definition of folders
folderin = fullfile(session.input_path,"Processed_DATA",ManipName);
folderout = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Matching','Rays');

%% Creation of folder containing rays data for parallele matching
if ~isfolder(folderout)
    mkdir(folderout);
end

% Rays data loading
fprintf("Rays loading...")
load(fullfile(folderin,'rays.mat'),'datacam')

Size = size(datacam(1).data);
nframes= Size(2);

% write results in file
fprintf("Writing data in process...\n")
for kframe=1:nframes
    if rem(kframe-1,nbFramePerJobMatching)==0
        if kframe~=1
            fclose(fid);
        end
        fileRays = fullfile(folderout,['rays_' num2str(kframe) '-' num2str(kframe+nbFramePerJobMatching-1) '.dat']);
        fprintf("%s\n",fileRays);
        fid=fopen(fileRays,'w');
    end
    Nrays=0;
    for kcam=1:numel(camID)
        Nrays = Nrays + numel(datacam(kcam).data(kframe).rayID);
    end
    fwrite(fid,Nrays,'uint32');
    
    for kcam=1:numel(camID)
        for kray=1:numel(datacam(kcam).data(kframe).rayID)
            fwrite(fid,camID(kcam),'uint8');
            fwrite(fid,datacam(kcam).data(kframe).rayID(kray),'uint16');
            fwrite(fid,datacam(kcam).data(kframe).P(kray,:),'float32');
            fwrite(fid,datacam(kcam).data(kframe).V(kray,:),'float32');
        end
    end
end
fclose(fid);
