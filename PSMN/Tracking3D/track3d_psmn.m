function [tracks,traj]=track3d_psmn(session,ManipName,NbFramePerJobMatching,minframe,maxframe,maxdist,lmin,flag_pred,npriormax,flag_conf)

% Load matches files created by parallelized matching jobs, make a single
% data array and call the function track3d_manualfit to track particles.
% --------------------------------------------------------------------------
% INPUTS
% session                : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
% ManipName              : Name of the experiment
% NbFramePerJobMatching  : Number of frame per job used for parallelized matching
% minframe               : number of the first frame to treat
% maxframe               : number of the last frame to treat. Pay attention, min and max
% frame have to be multiple of NbFramePerJob
% maxdist                : maximum travelled distance between two successive frames
% lmin : minimum length of a trajectory (number of frames)
% flag_pred              : 1 for predictive tracking, 0 otherwise
% npriormax              : maximum number of prior frames used for predictive tracking
% flag_conf              : 1 for conflict solving, 0 otherwise
%
% OUTPUTS
% traj(kt).ntraj : trajectory index
% traj(kt).L : trajectory length
% traj(kt).frames : trajectory frames
% traj(kt).x : x-position
% traj(kt).y : y-position
% traj(kt).z : z-position
% traj(kt).nmatch : element indices in tracks
% tracks : trajectory raw data
% --------------------------------------------------------------------------
% 2020-2021 D. Dumont
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if (maxframe-minframe+1) is a multiple of NbFramePerJob
if rem((maxframe-minframe+1),NbFramePerJobMatching)~=0
fprintf("(maxframe-minframe+1) has to be a multiple of NbFramePerJobTracking! Please correct it.")
return
end

%% Folder creation
folderin = fullfile(session.input_path,'Processed_DATA',ManipName,'Parallel','Matching','Rays')
folderout = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Tracking','Tracks')
if ~isfolder(folderout)
    mkdir(folderout);
end

%% Creation of FileName
FileName = sprintf('%d-%d',minframe,maxframe);

fprintf('Loading matches...\n')
tic

%load matches We load all matched files which contains matched point
%between min and maxframe
for i=minframe:NbFramePerJobMatching:maxframe
    fprintf(" frame %d",i)
    filematches=fullfile(folderin,['rays_' num2str(i) '-' num2str(i+NbFramePerJobMatching-1) '_out_cpp.h5']);
    fprintf('Input file: %s\n',filematches)

    if i==minframe
        for j=1:NbFramePerJobMatching
            Dataset=sprintf('/frame%d_xyze',j);
            if j==1
                 dataset = h5read(filematches,Dataset);
                 data=cat(2,ones(length(dataset),1)*(j+minframe-1),dataset);
            else

                dataset = h5read(filematches,Dataset);
                data1=cat(2,ones(length(dataset),1)*(j+minframe-1),dataset);
                data=cat(1,data,data1);
            end
        end
    else
        for j=1:NbFramePerJobMatching
            Dataset=sprintf('/frame%d_xyze',j);
            dataset = h5read(filematches,Dataset);
            data1=cat(2,ones(length(dataset),1)*(j+i-1),dataset);
            data=cat(1,data,data1);
        end
    end
end
elapsed=toc;
disp(['Elapsed time for data loading: ',num2str(elapsed,'%.2f'),'s']);
%% Call for track3d_manualfit function
[tracks,traj] = track3d_manualfit(folderout,FileName,data,maxdist,lmin,flag_pred,npriormax,flag_conf);

% %% Call for track3d_polyfit function
% [tracks,traj] = track3d_polyfit(folderout,FileName,data,maxdist,lmin,flag_pred,npriormax,flag_conf);

end