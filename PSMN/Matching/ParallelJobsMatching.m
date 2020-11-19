function ParallelJobsMatching(session,STM_path,ManipName,nframes,NbFramePerJobMatching,CamMatch,MaxDistance,nx,ny,nz,MaxMatchesPerRay,bminx,bmaxx,bminy,bmaxy,bminz,bmaxz,MinDistMatchperRay,Queue)
%% Create a file <ManipName>-ParallelMatching.sh which allows us to run all matching jobs on the PSMN with only one command.
% Create a folder Parallel with two subfolders SH and LOG to store sh and
% log file for each job.
%----------------------------------------------------------------------------
%%% Parameters : 
%%%     session                : Path to the achitecture root
%%%     ManipName              : Name of the folder experiment
%%%     nframes                : Total number of frames in the experiment
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test if Queue is defined or not
if ~exist('Queue','var')
   Queue = 'PIV'; 
end

%% Definition of folders
folderin = sprintf("%sProcessed_DATA/%s",session.input_path,ManipName);
folderout = sprintf("%sProcessed_DATA/%s",session.output_path,ManipName);
folderout_LOG = sprintf('%sProcessed_DATA/%s/Parallel/Matching/LOG',session.output_path,ManipName);
folderout_SH = sprintf('%sProcessed_DATA/%s/Parallel/Matching/SH',session.output_path,ManipName);
folderSTM = sprintf(STM_path);
RaysFolder = sprintf("%s/Parallel/Matching/Rays",folderin);

% Creation of folder containing .log files for parallel matching
if ~isfolder(folderout_LOG)
    mkdir(folderout_LOG);
end

% Creation of folder containing .sh files for parallel matching
if ~isfolder(folderout_SH)
    mkdir(folderout_SH);
end

for kframe=1:nframes
    if rem(kframe-1,NbFramePerJobMatching)==0
        fid = fopen(sprintf("%s/rays_%d-%d.sh",folderout_SH,kframe,kframe+NbFramePerJobMatching-1),'w');
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
        if Queue=='PIV'
            fwrite(fid,sprintf("#$ -q piv_debian*\n"));
            fwrite(fid,sprintf("#$ -P PIV\n"));
        else
            fwrite(fid,sprintf("#$ -q %s\n",Queue));
        end
        
        % Load user's environement for SGE
        fwrite(fid,sprintf("#$ -cwd\n"));
        % Exportation of environement variables over all execution nodes
        fwrite(fid,sprintf("#$ -V\n"));

        % Go into the work/submission directory. Important because without
        % it, the job is run from ~/.
        fwrite(fid,sprintf("cd ${SGE_O_WORKDIR}\n"));

        % Job execution
         fwrite(fid,sprintf("%sSTM -i %s/rays_%d-%d.dat -o %s/Parallel/Matching/Rays/ -f %d -c %d -d %f -m %d -x %d -y %d -z %d -b %d %d %d %d %d %d -s %d --hdf5> %s/rays_%d-%d.log\n",...
             folderSTM,RaysFolder,kframe, kframe+NbFramePerJobMatching-1,folderout, nframes,CamMatch, MaxDistance,MaxMatchesPerRay, nx,ny,nz,bminx,bmaxx,bminy,bmaxy,bminz,bmaxz, MinDistMatchperRay, folderout_LOG,kframe,kframe+NbFramePerJobMatching-1));

        fclose(fid);

    end
end

% Creation of file running all jobs: <ManipName>-paralleleMatching.sh
fid1=fopen(sprintf("%s/Parallel/%s-ParallelMatching.sh",folderout,ManipName),'w');
for kframe=1:nframes
    if rem(kframe-1,NbFramePerJobMatching)==0
        fwrite(fid1,sprintf("qsub %s/rays_%d-%d.sh\n",folderout_SH,kframe,kframe+NbFramePerJobMatching-1));
    end
end

fclose(fid1);