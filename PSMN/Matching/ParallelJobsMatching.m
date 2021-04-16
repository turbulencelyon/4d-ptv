function ParallelJobsMatching(session,STM_path,ManipName,minFrame,maxFrame,NbFramePerJobMatching,CamMatch,MaxDistance,nx,ny,nz,MaxMatchesPerRay,bminx,bmaxx,bminy,bmaxy,bminz,bmaxz,MinDistMatchperRay,Queue)
%% Create a file <ManipName>-ParallelMatching.sh which allows us to run all matching jobs on the PSMN with only one command.
% Create a folder Parallel with two subfolders SH and LOG to store sh and
% log file for each job.
%----------------------------------------------------------------------------
%%% Parameters : 
%%%     session                : Path to the achitecture root
%%%     STM_path               : Path to the STM script
%%%     ManipName              : Name of the folder experiment
%%%     minFrame               : number of the first frame in the experiment
%%%     maxFrame               : number of the last frame in the experiment
%%%     NbFramePerJob          : Number of frames per job
%%%     CamMatch               : Minimum number of rays to get a match
%%%     MaxDistance            : Maximal authorized distance between rays to
%%%     consider having a match
%%%     nx,ny,nz               : number of voxels in each direction
%%%     MaxMatchesPerRay       : Maximum number of matches for one ray. 2 to
%%%     consider particle overlap
%%%     NbFramePerJobMatching  : number of frames per job. Pay attention, has to be
%%%     chosen as a function of processing time of one picture, in order to
%%%     that each job runs for 10 min (PSMN requirements).
%%%     bminx=min distance in x
%%%     bmaxx=max distance in x 
%%%     bminy=min distance in y
%%%     bmaxy=max distance in y 
%%%     bminz=min distance in z
%%%     bmaxz=max distance in z 
%%%     MaxDistMatchperRay= specifie a volume in wich you cannot have an
%%%     other match if you already found one (avoid to consider several
%%%     match for one particule)
%%%     Queue      s           : Running queue. By defaut it is equal to 'PIV'.
%%%     It is possible to run jobs on monointeldeb128 or monointeldeb48 for
%%%     example. Do 'qstat -g c' to get all opened queues.
%------------------------------------------------------------------------------
% 2020-2021 D. Dumont
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test if Queue is defined or not
if ~exist('Queue','var')
   Queue = 'PIV'; 
end

%% Definition of folders
folderin = fullfile(session.input_path,"Processed_DATA",ManipName);
folderout = fullfile(session.output_path,"Processed_DATA",ManipName);
folderout_LOG = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Matching','LOG');
folderout_SH = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Matching','SH');
folderSTM = fullfile(STM_path);
RaysFolder = fullfile(folderin,"Parallel","Matching","Rays");

% Creation of folder containing .log files for parallel matching
if ~isfolder(folderout_LOG)
    mkdir(folderout_LOG);
end

% Creation of folder containing .sh files for parallel matching
if ~isfolder(folderout_SH)
    mkdir(folderout_SH);
end

% Creation of folder containing .log files for parallel tracking
folderout_tracking_LOG = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Tracking','LOG');
folderout_tracking_OUT = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Tracking','OUT');
if ~isfolder(folderout_tracking_LOG)
    mkdir(folderout_tracking_LOG);
end

% Creation of folder containing .sh files for parallel tracking
if ~isfolder(folderout_tracking_OUT)
    mkdir(folderout_tracking_OUT);
end

% Creation of folder containing .log files for parallel stitching
folderout_stitching_LOG = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Stitching','LOG');
folderout_stitching_OUT = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Stitching','OUT');
if ~isfolder(folderout_stitching_LOG)
    mkdir(folderout_stitching_LOG);
end

% Creation of folder containing .sh files for parallel stitching
if ~isfolder(folderout_stitching_OUT)
    mkdir(folderout_stitching_OUT);
end

for kframe=minFrame:NbFramePerJobMatching:maxFrame-NbFramePerJobMatching+1
    if rem(kframe-minFrame,NbFramePerJobMatching)==0
        fid = fopen(fullfile(folderout_SH,['rays_' num2str(kframe) '-' num2str(kframe+NbFramePerJobMatching-1) '.sh']),'w');
        fwrite(fid,sprintf("#!/bin/bash\n"));
        
        %%% SGE variables
        %%% job shell
        fwrite(fid,sprintf("#$ -S /bin/bash\n"));
        
        %%% job name
        fwrite(fid,sprintf("#$ -N rays_%d-%d\n",kframe, kframe+NbFramePerJobMatching-1));
        fwrite(fid,sprintf("#$ -e %s/rays_%d-%d.log\n",folderout_LOG,kframe, kframe+NbFramePerJobMatching-1));
        fwrite(fid,sprintf("#$ -o %s/rays_%d-%d.log\n",folderout_LOG,kframe, kframe+NbFramePerJobMatching-1));

        %%% Waiting queue
%         fwrite(fid,sprintf("#$ -q piv_debian* ## !(*test*|*gpu*)\n"));
        if strcmp(Queue,'PIV')
            fwrite(fid,sprintf("#$ -q piv_debian*\n"));
            fwrite(fid,sprintf("#$ -P PIV\n"));
        else
            fwrite(fid,sprintf("#$ -q %s\n",Queue));
        end
        
        % Load user's environment for SGE
        fwrite(fid,sprintf("#$ -cwd\n"));
        % Exportation of environement variables over all execution nodes
        fwrite(fid,sprintf("#$ -V\n"));

        % Go into the work/submission directory. Important because without
        % it, the job is run from ~/.
        fwrite(fid,sprintf("cd ${SGE_O_WORKDIR}\n"));

        % Job execution
         fwrite(fid,sprintf("%sSTM -i %s/rays_%d-%d.dat -o %s/Parallel/Matching/Rays/ -f %d -c %d -d %f -m %d -x %d -y %d -z %d -b %d %d %d %d %d %d -s %d --hdf5> %s/rays_%d-%d.log\n",...
             folderSTM,RaysFolder,kframe, kframe+NbFramePerJobMatching-1,folderout, NbFramePerJobMatching, CamMatch, MaxDistance,MaxMatchesPerRay, nx,ny,nz,bminx,bmaxx,bminy,bmaxy,bminz,bmaxz, MinDistMatchperRay, folderout_LOG,kframe,kframe+NbFramePerJobMatching-1));

        fclose(fid);

    end
end
if kframe+NbFramePerJobMatching-1~=maxFrame
    fid = fopen(fullfile(folderout_SH,['rays_' num2str(kframe+NbFramePerJobMatching) '-' num2str(maxFrame) '.sh']),'w');
    fwrite(fid,sprintf("#!/bin/bash\n"));

    %%% SGE variables
    %%% job shell
    fwrite(fid,sprintf("#$ -S /bin/bash\n"));

    %%% job name
    fwrite(fid,sprintf("#$ -N rays_%d-%d\n", kframe+NbFramePerJobMatching, maxFrame));
    fwrite(fid,sprintf("#$ -e %s/rays_%d-%d.log\n",folderout_LOG,kframe+NbFramePerJobMatching, maxFrame));
    fwrite(fid,sprintf("#$ -o %s/rays_%d-%d.log\n",folderout_LOG,kframe+NbFramePerJobMatching, maxFrame));

    %%% Waiting queue
%         fwrite(fid,sprintf("#$ -q piv_debian* ## !(*test*|*gpu*)\n"));
    if strcmp(Queue,'PIV')
        fwrite(fid,sprintf("#$ -q piv_debian*\n"));
        fwrite(fid,sprintf("#$ -P PIV\n"));
    else
        fwrite(fid,sprintf("#$ -q %s\n",Queue));
    end

    % Load user's environment for SGE
    fwrite(fid,sprintf("#$ -cwd\n"));
    % Exportation of environement variables over all execution nodes
    fwrite(fid,sprintf("#$ -V\n"));

    % Go into the work/submission directory. Important because without
    % it, the job is run from ~/.
    fwrite(fid,sprintf("cd ${SGE_O_WORKDIR}\n"));

    % Job execution
     fwrite(fid,sprintf("%sSTM -i %s/rays_%d-%d.dat -o %s/Parallel/Matching/Rays/ -f %d -c %d -d %f -m %d -x %d -y %d -z %d -b %d %d %d %d %d %d -s %d --hdf5> %s/rays_%d-%d.log\n",...
         folderSTM,RaysFolder,kframe+NbFramePerJobMatching,maxFrame,folderout, NbFramePerJobMatching, CamMatch, MaxDistance,MaxMatchesPerRay, nx,ny,nz,bminx,bmaxx,bminy,bmaxy,bminz,bmaxz, MinDistMatchperRay, folderout_LOG,kframe+NbFramePerJobMatching,maxFrame));

    fclose(fid);
end

% Creation of file running all jobs: <ManipName>-paralleleMatching.sh
fid1=fopen(fullfile(folderout,'Parallel',[char(ManipName) '-ParallelMatching.sh']),'w');
for kframe=minFrame:maxFrame
    if rem(kframe-1,NbFramePerJobMatching)==0
        fwrite(fid1,sprintf("qsub %s/rays_%d-%d.sh\n",folderout_SH,kframe,kframe+NbFramePerJobMatching-1));
    end
end
if kframe+NbFramePerJobMatching-1~=maxFrame
    fwrite(fid1,sprintf("qsub %s/rays_%d-%d.sh\n",folderout_SH,kframe+NbFramePerJobMatching,maxFrame));
end
    
fclose(fid1);