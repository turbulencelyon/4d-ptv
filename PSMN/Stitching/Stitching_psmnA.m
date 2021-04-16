function StitchedTraj = Stitching_psmnA(session,ManipName,minframe,maxframe,dfmax,dxmax,dvmax,lmin,Test)
% Create automatic filepath, load data and call for stitchTracks function.
% To use after track3d_psmn.m function, only for the first run.
% Reconnect trajectories in the whole data.
%----------------------------------------------------------------------------------------
% Parameters:
%   session                : session.path contains MyPath,
%   ManipName              : Name of the experiment,
%   minframe               : First frame to process,
%   maxframe               : Last frame to process,
%   dfmax                  : maximum number of tolerated missing frames to
%   reconnect to trajectories,
%   dxmax                  : maximum tolerated distance (in norm) between projected
%   point after the first trajectory and the real beginning position of the
%   stitched one,
%   dvmax                  : maximum tolerated relative velocity difference between
%   the end of the first trajectory and the beginning of the stitched one,
%   lmin                   : minimum length for a trajectory to be stitched.
% ----------------------------------------------------------------------------------------
% 04/2020 - David Dumont
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic()
%% Input and output folders and filepath
folderin = fullfile(session.input_path,'Processed_DATA', ManipName, 'Parallel','Tracking','Tracks')
folderout = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Stitching','StitchA')
if exist('Test','var')
    filepath = fullfile(folderout,['StitchedTracksA_' num2str(minframe) '-' num2str(maxframe) '_dfmax' num2str(dfmax) '_' num2str(Test)]);
else
    filepath = fullfile(folderout,['StitchedTracksA_' num2str(minframe) '-' num2str(maxframe) '_dfmax' num2str(dfmax)]);
end

%% Creation of output folder if it does not exist
if ~isfolder(folderout)
    mkdir(folderout)
end  

fprintf("Data loading")
tic()
FileName = fullfile(folderin,['tracks_' num2str(minframe) '-' num2str(maxframe)])
traj = h52tracks(FileName);
toc()

%% Let's stitch trajectories everywhere
StitchedTraj = stitchTracks(traj,dfmax,filepath,dxmax,dvmax,lmin);
toc()