function [tracks,traj]=track3d(session,ManipName,FileName,maxdist,lmin,flag_pred,npriormax,flag_conf)

% Load matches data and give it to track3d_manualfit to track particles.
% ____________________________________________________________________________
% INPUTS
% session   : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
% ManipName : Name of the experiment
% FileName  : Name of the matched file without its extension (without .h5)
% maxdist   : maximum travelled distance between two successive frames
% lmin      : minimum length of a trajectory (number of frames)
% flag_pred : 1 for predictive tracking, 0 otherwise
% npriormax : maximum number of prior frames used for predictive tracking
% flag_conf : 1 for conflict solving, 0 otherwise
%
% OUTPUTS
% traj(kt).ntraj  : trajectory index
% traj(kt).L      : trajectory length
% traj(kt).frames : trajectory frames
% traj(kt).x      : x-position
% traj(kt).y      : y-position
% traj(kt).z      : z-position
% traj(kt).nmatch : element indices in tracks
% tracks          : trajectory raw data
% ____________________________________________________________________________

% Folder name creation
folderin = sprintf("%sProcessed_DATA/%s/",session.input_path,ManipName);
folderout = sprintf("%sProcessed_DATA/%s/",session.output_path,ManipName);

if ~isfolder(folderout)
    makedirs(folderout)
end


disp('Loading matches...');
tic
data = h5read(sprintf("%s%s.h5",folderin,FileName),'/tracks');

%% Call for track3d_manualfit function
[tracks,traj] = track3d_manualfit(folderout,FileName,data,maxdist,lmin,flag_pred,npriormax,flag_conf);

% %% Call for track3d_polyfit function
% [tracks,traj] = track3d_polyfit(folderout,FileName,data,maxdist,lmin,flag_pred,npriormax,flag_conf);
end