function [P,V]=Centers2Rays(session,ManipName,Calib,camID,Ttype)
%%   To DO     - h5 output
%%%   For a given set of calibration points and center files with
%%%   pixel positions, determines the corresponding rays of light.
%%%   Save rays into a .mat file and in rays.dat file.
%----------------------------------------------------------------------------------------
%%% INPUT Parameters:
%%%   session     : Path to the achitecture root (2 fields: session.input_path
%%% and session.output_path)
%%%   ManipName   : Name of the folder experiment
%   Calib       : Path of the calibration file
%   camID       : List of camera numbers. ex: [1,2,3] if you have 3 cameras numbered 1,2,3 respectively.
%   Ttype       : Type of the calibration used. 'T1'-> Linear
%   transformation, 'T3' -> Cubic transformation.
% ------------------------------------------------------------------------------------------
% 2020-2021 : D. Dumont (adapted from M. Bourgoin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If Ttype does not exist it is set to T1 (linear transformation)
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
fprintf("Loading calibration file... \n")
load(Calib,'calib');

% Open results ray file
fileRays=fullfile(folderout,'rays.dat');

fid=fopen(fileRays,'w');

fprintf("Let's start loop over cameras...\n")
tic()
% Loop over cameras
for kcam=1:numel(camID)
    % select calibration for camera camID(kcam)

    % new version without interpolant
    calibNcam=calib(:,kcam);
    
    % load centers for camera camID(kcam)
    fileCenters=fullfile(folderin,['centers_cam' num2str(camID(kcam)) '.mat']);
    [CC,firstFrame,endFrame] = readCentersMAT(fileCenters); 
    
    % loop over frames
    for k=firstFrame:endFrame
        if rem(k,100)==0
            fprintf("cam %d frame %d...\n",kcam,k)
        end
        % convert pixel coordinates into rays of light using the
        % calibration
        if ~isempty(CC(k).X)
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
        else
                data(k).P=[];
                data(k).V=[];
                data(k).rayID=[];
        end
    end
    datacam(kcam).data=data;
end
toc()

% Writing data in matlab file
fprintf("Writing data in matlab file in progress...\n")
save(fullfile(folderout,'rays.mat'),'datacam','Calib','-v7.3')

%%
% write results in file
fprintf("Writing data in process...\n")
for kframe=firstFrame:endFrame
    Nrays=0;
    for kcam=1:numel(camID)
        Nrays = Nrays + numel(datacam(kcam).data(kframe).rayID);
    end
    fwrite(fid,Nrays,'uint32');

    for kcam=1:numel(camID)
        for kray=1:numel(datacam(kcam).data(kframe).rayID)
            fwrite(fid,camID(kcam),'uint8');
            fwrite(fid,datacam(kcam).data(kframe).rayID(kray),'uint16');
            fwrite(fid,datacam(kcam).data(kframe).P(kray,:),'float32');
            fwrite(fid,datacam(kcam).data(kframe).V(kray,:),'float32');
        end
    end
end

fclose(fid);
