function CalibTest(Calib_path,nb_plan)

%% This function draw the rays between the calibration point and the camera
% Input:
% Calib_path        : the path of the calibration file calib.mat
% nb_plan           : the number of calibration plane you used 
%
% Ouput:
% A file rays_test_particle.mat 


% If Ttype does not exist it is set to T1 (linear transformation)
if ~exist('Ttype','var')
    Ttype='T1';
end

A=open(Calib_path);
calib=A.calib;

for kcam=1:4
    % select calibration for camera camID(kcam)
    % old version with interpolant
    % calibNcam=calibInterp(camID(kcam)+1).calibInterp;
    
    % new version without interpolant
    calibNcam=calib(:,kcam);
    
    % load centers for the calibration
    Xtt=[];
    Ytt=[];
    for k=1:nb_plan
%         X= calib(k,kcam).pimg(:,1); 
%         Y= calib(k,kcam).pimg(:,2); 
        fileCenters=sprintf("/Xnfs/convection/Stage_EB_2020/Calibration_Test/center_cam%d.mat", kcam);
        CC = readCentersMAT(fileCenters); 
    % Compute the rays 
  

        % convert pixel coordinates into rays of light using the
        % calibration
        [P,V]=findRays(calibNcam,CC.X',CC.Y',Ttype);
        
        % exclude particles for which rays are obtained by extrapolation
        % outside the actually calibrated convex hull
        rayID=find((~isnan(P(:,1)))&(~isempty(P(:,1))));
        if ~isempty(rayID)
            data(k).P=P(rayID,:);
            data(k).V=V(rayID,:);
            data(k).rayID=rayID;
        end
    datacam(kcam).data=data;
    end
end
save(sprintf('/Xnfs/convection/Stage_EB_2020/Processed_DATA/Ra1.60e10_peudense_2/rays_test_particules.mat'),'datacam','-v7.3')
end
