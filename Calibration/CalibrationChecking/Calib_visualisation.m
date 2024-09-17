function Calib_visualisation(dirIn,filemane,CamID,Nplane)

%% Plot the calibration points on the raw image in order to check if the calibration is good 
%% Input 
% dirIn     : the directory where the raw image and the calib file are
% Nplane    : the number of the plane you want to check
%% Output 
% N figure (where N is the number of cameras) where the cross are the point
% you place during the claibration
%The idea is to check that the point are centered on the calibration plate
% /!\ It's really usefull to check because bad calibration=bad result
%% 


%%Load the calibration and definition of pimg
calib_path=fullfile(dirIn,filemane);
load(calib_path);


%%Load the image file
for kcam=CamID
    PimgX=calib(Nplane,kcam).pimg(:,1);
    PimgY=calib(Nplane,kcam).pimg(:,2);
    figure('numberTitle','off','Name',sprintf('Cam %d',kcam))
    filename = sprintf('%s/CalibrationPlan_%d_cam%d.%s',dirIn, Nplane, kcam,'tif');        
    Img=imread(filename);
    imshow(Img);
    hold on 
    plot(PimgX,PimgY,'rx',LineWidth=3);
    PimgX=[];
    PimgY=[];
    filename=[];
end
end 
