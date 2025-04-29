function [tracks,traj]=track3d(session,ManipName,FileName,NbFrame,maxdist,lmin,flag_pred,npriormax,flag_conf,minFrame,test)
%function [tracks,traj]=track3d(session,ManipName,FileName,nCam,NbFrame,maxdist,lmin,flag_pred,npriormax,flag_conf,minFrame,test)

% Load matches data and give it to track3d_manualfit to track particles.
% ____________________________________________________________________________
% INPUTS
% session   : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
% ManipName : Name of the experiment
% FileName  : Name of the matched file without its extension (without .h5)
%%% nCam      : Number of camera in the used for rays calculation and matches
% NbFrame   : Number of frame in the file
% maxdist   : maximum travelled distance between two successive frames
% lmin      : minimum length of a trajectory (number of frames)
% flag_pred : 1 for predictive tracking, 0 otherwise
% npriormax : maximum number of prior frames used for predictive tracking
% flag_conf : 1 for conflict solving, 0 otherwise
% minFrame : (optional) number of the first frame. Default = 1.
% test      : (optional) allows you to not save data when you are doing
% tests to find best parameters
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
%
% 2020-2021 D. Dumont (adapted from B. Bourgoin)
% ____________________________________________________________________________

if ~exist('minFrame','var')
    minFrame=1;
end
if ~exist('test','var')
    test=false;
end

% Folder name creation
folderin = fullfile(session.input_path, 'Processed_DATA', ManipName);
folderout = fullfile(session.output_path, 'Processed_DATA', ManipName);

if ~isfolder(folderout)
    mkdir(folderout)
end

% disp('Loading matches...');
filename = fullfile(folderin,FileName);

data = h52matchesTracking(filename,NbFrame,minFrame);
%data = h52matches(filename,NbFrame,minFrame,nCam);
%% Call for track3d_manualfit function
[tracks,traj] = track3d_manualfit(folderout,FileName,data,maxdist,lmin,flag_pred,npriormax,flag_conf,test);

% %% Call for track3d_polyfit function
% [tracks,traj] = track3d_polyfit(folderout,FileName,data,maxdist,lmin,flag_pred,npriormax,flag_conf);
end
