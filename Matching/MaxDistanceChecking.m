function MaxDistanceChecking(session,ManipName,nframe,FirstFrame,nbins,FilePath)
% Plot the matching error histogram to check if the maxdistance value is
% correct or not
% -------------------------------------------------------------
% INPUT
% session               : Path to the achitecture root (2 fields: session.input_path and session.output_path)
% ManipName             : Name of the folder experiment
% nframe                : number of frames treated during the test
% FirstFrame (optional) : number of the first frame (=1 by defaut)
% nbins (optional)      : number of histogram bins (=100 by defaut)
% FilePath (optional)   : path of the matches file
% -------------------------------------------------------------
% OUTPUT
% Matching error histogram
% -------------------------------------------------------------
% 2021 D. Dumont
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if FirstFrame does not exist we create it
if ~exist('FirstFrame','var')
    FirstFrame = 1;
end

%% if nbins does not exist we create it
if ~exist('nbins','var')
    nbins = 100;
end

%% Construction of FilePath if it does not exist
if ~exist('FilePath','var')
    FilePath = fullfile(session.input_path,'Processed_DATA',ManipName,'Parallel','Matching','Rays',['rays_' num2str(FirstFrame) '-' num2str(nframe) '_out_cpp'])
end

%% Loading of matches
matches = h52matches(FilePath,nframe,FirstFrame);

%% Plot of histogram
figure
histogram(matches(:,5),nbins)
xlabel('Matching error')
ylabel('Histogram')