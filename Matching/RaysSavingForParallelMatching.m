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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Definition of folders
folderin = sprintf("%s/Processed_DATA/%s",session.input_path,ManipName);
folderout = sprintf('%s/Processed_DATA/%s/Parallel/Matching/Rays',session.output_path,ManipName);

%% Creation of folder containing rays data for parallele matching
if ~isfolder(folderout)
    mkdir(folderout);
end

% Rays data loading
fprintf("Rays loading...")
load(sprintf('%s/rays.mat',folderin),'datacam')

Size = size(datacam(1).data);
nframes= Size(2);

% write results in file
fprintf("Writing data in process...\n")
for kframe=1:nframes
    if rem(kframe-1,nbFramePerJobMatching)==0
        if kframe~=1
            fclose(fid);
        end
        fileRays = sprintf("%s/rays_%d-%d.dat",folderout,kframe,kframe+nbFramePerJobMatching-1);
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
