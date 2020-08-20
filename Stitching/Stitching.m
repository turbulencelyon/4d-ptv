function StitchedTraj = Stitching(session,ManipName,FileName,dfmax,dxmax,dvmax,lmin,window_length,poly_order,framerate)
% Create automatic filepath, load data and call for stitchTracks function.
% To use after track3d.m function.
% 04/2020 - David Dumont
%----------------------------------------------------------------------------------------
% Parameters:
%   session        : session.path contains MyPath, (2 fields: session.input_path
% and session.output_path),
%   ManipName      : Name of the experiment,
%   FileName       : Name of the tracks file,
%   dfmax          : maximum number of tolerated missing frames to
%   reconnect to trajectories,
%   dxmax          : maximum tolerated distance (in norm) between projected
%   point after the first trajectory and the real beginning position of the
%   stitched one,
%   dvmax          : maximum tolerated relative velocity difference between
%   the end of the first trajectory and the beginning of the stitched one,
%   lmin           : minimum length for a trajectory to be stitched,
% window_length : size of the window used for savgol filtering (to compute speed),
% poly_order    : order of the polynom used for savgol filtering (to
% compute velocity),
% framerate     : framerate of the experiment.
% ------------------------------------------------------------------------------------------

datapath = sprintf('%sProcessed_DATA/%s/%s',session.input_path,ManipName,FileName);
filepath = sprintf('%sProcessed_DATA/%s/stiched_%s_dfmax%d',session.output_path,ManipName,FileName,dfmax);

%% Data loading
traj = h52tracks(datapath);

fprintf("Before stitching, we have %d trajectories\n",numel(traj))

%% Let's stitch trajectories !
StitchedTraj = stitchTracks(traj,dfmax,filepath,dxmax,dvmax,lmin);

end

