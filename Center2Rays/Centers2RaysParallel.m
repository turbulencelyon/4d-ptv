function [P,V]=Centers2RaysParallel(session,ManipName,Calib,kcam,FirstFrame,Ttype)
%%% The same same function than Centers2Rays but for one camera (allow to run the n camera in paralell)
%%%       
%%%   For a given set of calibration points and center files with
%%%   pixel positions, determines the corresponding rays of light.
%%%----------------------------------------------------------------------------------------
%%% INPUT Parameters:
%%%   session     : Path to the achitecture root (2 fields: session.input_path
%%% and session.output_path)
%%%   ManipName   : Name of the folder experiment
%%%   Calib       : Path of the calibration file
%%%   kcam        : The number of the camera one which you want to compute the
%%%   rays 
%%%   FirstFrame  : The first frame you want to treat 
%%%   Ttype       : Type of the calibration used. 'T1'-> Linear
%%%   transformation, 'T3' -> Cubic transformation.
%%% ------------------------------------------------------------------------------------------
%%% OUTPUT Parameters:
%%% A file 'rays_camX.mat'. Once you have the rays for all the camera, use
%%% the function 'Rays_recombinaison.m' to create the complete file rays.mat 
%%% ------------------------------------------------------------------------------------------
%%% To compile use the command:
%%% mcc -m submission_Center2Rays.m -a /applis/PSMN/generic/Matlab/R2019b/toolbox/images/images
%%%
%%% If Ttype does not exist it is set to T1 (linear transformation)
%%% ------------------------------------------------------------------------------------------
% 2020-2021 : E. Bernard (adapted from M. Bourgoin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('Ttype','var')
    Ttype='T1';
end

%% Definition of folders
folderin = fullfile(session.input_path, 'Processed_DATA', ManipName);
folderout = fullfile(session.output_path, 'Processed_DATA', ManipName);
if ~isfolder(folderout)
    mkdir(folderout);
end

% Load calibration file (which contains calibration transformations for all
% cameras).
%% Charger la calib
fprintf("Loading calibration file %s... \n",Calib)
load(Calib,'calib');


fprintf("Calibration loading...\n")
% select calibration for camera camID(kcam)
% new version without interpolant
calibNcam=calib(:,kcam);


% load centers for camera camID(kcam)
fileCenters=fullfile(folderin,['centers_cam' num2str(kcam) '.mat']);
[CC,nframes] = readCentersMAT(fileCenters); 

fprintf("Let's loop over all frames\n")
% loop over frames
for k=FirstFrame:nframes

    if rem(k,100)==0
        fprintf("cam %d frame %d...\n",kcam,k)
    end
    % convert pixel coordinates into rays of light using the
    % calibration
    [P,V]=findRays(calibNcam,CC(k).X',CC(k).Y',Ttype);

    % exclude particles for which rays are obtained by extrapolation
    % outside the actually calibrated convex hull
    rayID=find((~isnan(P(:,1)))&(~isempty(P(:,1))));
    if ~isempty(rayID)
        data(k).P=P(rayID,:);
        data(k).V=V(rayID,:);
        data(k).rayID=rayID;
%             kframes = kframes + 1;
    end
end
datacam(kcam).data=data;



% Writing data in matlab file
fprintf("Writing data in matlab file in progress...\n")
save(fullfile(folderout,['rays_cam' num2str(kcam) '.mat']),'datacam','-v7.3')

%% The results have then to be merged manualy with the comment code above
