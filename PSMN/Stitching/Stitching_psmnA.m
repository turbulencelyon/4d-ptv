function StitchedTraj = Stitching_psmnA(session,ManipName,minframe,maxframe,dfmax,dxmax,dvmax,lmin)
% Create automatic filepath, load data and call for stitchTracks function.
% To use after track3d_psmn.m function, only for the first run.
% Reconnect trajectories in the whole data.
% 04/2020 - David Dumont
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
% ------------------------------------------------------------------------------------------
tic
%% Input and output folders and filepath
folderin = fullfile(session.input_path,'Processed_DATA', ManipName, 'Parallel','Tracking','Tracks');
folderout = fullfile(session.output_path,'Processed_DATA',ManipName,'Parallel','Stitching','StitchA');
filepath = fullfile(folderout,['StitchedTracksA_' num2str(minframe) '-' num2str(maxframe) '_dfmax' num2str(dfmax)]);

%% Creation of output folder if it does not exist
if ~isfolder(folderout)
    mkdir(folderout)
end  

traj = h52tracks(fullfile(folderin,['tracks_' num2str(minframe) '-' num2str(maxframe)]));

%% Let's stitch trajectories everywhere
StitchedTraj = stitchTracks(traj,dfmax,filepath,dxmax,dvmax,lmin);
toc