function StitchedTraj = Stitching_psmnB(session,ManipName,minframe,maxframe,NbFramePerJobTracking,dfmax,dxmax,dvmax,lmin)
% Create automatic filepath, load data and call for stitchTracks function.
% To use after track3d_psmn.m function, only for the first run.
% Reconnect trajectories in the whole data.
% 04/2020 - David Dumont
%----------------------------------------------------------------------------------------
% Parameters:
%   session                : session.path contains MyPath, (2 fields: session.input_path
% and session.output_path),
%   ManipName              : Name of the experiment,
%   minframe               : First frame to process,
%   maxframe               : Last frame to process,
%   NbFramePerJobTracking  : Number of frame per tracking job,
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
folderin = sprintf('%sProcessed_DATA/%s/Parallel/Stitching/StitchA/',session.path,ManipName)
folderout = sprintf('%sProcessed_DATA/%s/Parallel/Stitching/',session.path,ManipName)
filepath = sprintf('%sStitchedTracksB_%d-%d_dfmax%d',folderout,minframe,maxframe,dfmax)


%% Creation of output folder if it does not exist
if ~isfolder(folderout)
    mkdir(folderout)
end  

%% Data loading
for i=minframe:NbFramePerJobTracking:maxframe
    if i==minframe
        traj = h52stitch(sprintf('%sStitchedTracksA_%d-%d_dfmax60',folderin,i,i+NbFramePerJobTracking-1))
    else
        Traj = h52stitch(sprintf('%sStitchedTracksA_%d-%d_dfmax60',folderin,i,i+NbFramePerJobTracking-1));
        traj_index = num2cell((numel(traj)+1):(numel(traj)+numel(Traj)));
        [Traj.ntraj] = traj_index{:};
        traj(end+1:end+numel(Traj)) = Traj;
    end
end

%% Let's stitch trajectories only at the files limits without saving output (filepath not given to stitchTracksSides)
StitchedTraj = stitchTracksSides(traj,minframe,maxframe,NbFramePerJobTracking,dfmax,dxmax,dvmax,lmin,filepath);
toc