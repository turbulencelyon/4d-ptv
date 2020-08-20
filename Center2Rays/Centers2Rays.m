function [P,V]=Centers2Rays(session,ManipName,Calib,camID,Ttype)
%% TO DO: - split this function in several to run as many jobs as camera and create an unique file later !!!
%         - h5 output
%   For a given set of calibration interpolanst and center files with
%   pixel positions, determines the corresponding rays of light.
%   camID is an array containing the IDs of the cameras center files to be
%   processed (for instance camID=[0 1 2] for a 3 cameras tracking
%   experiment).
%----------------------------------------------------------------------------------------
% INPUT Parameters:
%   session     : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
%   ManipName   : Name of the folder experiment
%   Calib       : Path of the calibration file
%   camID       : List of camera numbers. ex: [1,2,3] if you have 3 cameras numbered 1,2,3 respectively.
%   Ttype       : Type of the calibration used. 'T1'-> Linear
%   transformation, 'T3' -> Cubic transformation.
% ------------------------------------------------------------------------------------------

% If Ttype does not exist it is set to T1 (linear transformation)
if ~exist('Ttype','var')
    Ttype='T1';
end

%% Definition of folders
folderin = sprintf("%sProcessed_DATA/%s",session.input_path,ManipName);
folderout = sprintf("%sProcessed_DATA/%s",session.output_path,ManipName);

% Load calibration file (which contains calibration interpolants for all
% cameras).
%% Charger la calib
fprintf("Loading calibration file... \n")
load(Calib,'calib');

% Open results ray file
fileRays=sprintf('%s/rays.dat',folderout);
if ~isfolder(folderout)
    mkdir(folderout);
end

fid=fopen(fileRays,'w');

fprintf("Let's start loop over cameras...")
% Loop over cameras
for kcam=1:numel(camID)
    % select calibration for camera camID(kcam)
    % old version with interpolant
    % calibNcam=calibInterp(camID(kcam)+1).calibInterp;
    
    % new version without interpolant
    calibNcam=calib(:,kcam);
    
    % load centers for camera camID(kcam)
    fileCenters=sprintf("%s/centers_cam%d.mat",folderin, camID(kcam));
    [CC,nframes] = readCentersMAT(fileCenters); 
    
    % loop over frames
    for k=1:nframes
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
end

% Writing data in matlab file
fprintf("Writing data in matlab file in progress...")
save(sprintf('%s/rays.mat',folderout),'datacam','-v7.3')

%%
% write results in file
fprintf("Writing data in process...\n")
for kframe=1:nframes
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
