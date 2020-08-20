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
folderin = sprintf('%sProcessed_DATA/%s/Parallel/Tracking/Tracks/',session.input_path,ManipName);
folderout = sprintf('%sProcessed_DATA/%s/Parallel/Stitching/StitchA/',session.output_path,ManipName);
filepath = sprintf('%sStitchedTracksA_%d-%d_dfmax%d',folderout,minframe,maxframe,dfmax);

%% Creation of output folder if it does not exist
if ~isfolder(folderout)
    mkdir(folderout)
end  

traj = h52tracks(sprintf('%stracks_%d-%d',folderin,minframe,maxframe));

%% Let's stitch trajectories everywhere
StitchedTraj = stitchTracks(traj,dfmax,filepath,dxmax,dvmax,lmin);
toc