function [calib] = MakeCalibration(dirIn,zPlanes,camName,gridSpace,th,dotSize,lnoise,blackDots,extension,FirstPlane,FirstCam)
% Make calibration file.
% With human help, detect calibration grid on calibration pictures and
% compute 1st and 3rd order transformations.
% Store everything in a structure calib.mat saved in the same directory
%----------------------------------------------------------------------------------------
% INPUT Parameters:
%   dirIn                 : directory containning calibration pictures. calib.mat is saved in this directory,
%   zPlanes               : list of z position of each calibration plane. Given in SI units.
%   camName               : list of camera number if 3 cameras camName={1,2,3},
%   gridSpace             : distance between two points in the calibration mire. Given in SI units,
%   th                    : threshold to detect points on pictures,
%   dotSize               : typical diameter (in pixels) of points,
%   lnoise                : typical correlation size of noise (lnoise =1 by defaut) -> Use to make a gaussian filtering of pictures,
%   blackDots             : if true, script will look for black dots on white
%   background. If false script will look for white dots on dark background,
%   extension (optional)  : pictures extension. By defaut extention = 'tif', 
%   FirstPlane (optional) : number of the first plane to treat (useful
%   when you did a mistake during calibration to start at the right plane),
%   FirstCam (optional)   : number of the first camera to treat (idem).
% 
% OUTPUT
%     a calib.mat file saved in dirIn which contains the structure 'calib' such as for the kz plane and the kcam camera:
%     
%     calib(kz,kcam).posPlane        : z position which corresponds to kz plane
%     calib(kz,kcam).pimg            : detected mire points on the calibration picture (2D) in px units,
%     calib(kz,kcam).pos3D           : detected mire points in 3D space and in SI units,
%     calib(kz,kcam).movedpoints     : index of moved points,
%     calib(kz,kcam).addedpoints     : index of added points,
%     calib(kz,kcam).T1rx2px         : Linear transformation from real
%     world to px. T1rw2px=inverse(T1px2rw), Not saved as only one of the
%     two is sufficient.
%     calib(kz,kcam).T3rw2px         : Cubic transformation from real world to px. Pay attention T3rw2px~=inverse(T3px2rw) !,
%     calib(kz,kcam).T1px2rw         : Linear transformation from px to real world,
%     calib(kz,kcam).T3px2rw         : Cubic transformation from px to real world,
%     calib(kz,kcam).cHull           : coordinates of the region of interest,
%     calib(kz,kcam).name            : camera number (kcam),
%     calib(kz,kcam).dirPlane        : [i,j,k] axes orientation. i,j = axis provided by the calibration mire, k is the axis along which the mire is displaced. For instance, if the mire is in the plane (Oxy) and that it is displaced along z axis, dirPlane=[1,2,3]. 
    
% IMPORTANT! (the following is kind of confusing)
%
% To transform image to real world use:
% [x_rw,y_rw,z_rw] = transformPointsInverse(Trw2px,posimg(:,1),posimg(:,2));
%
% To transform real world positions to image positions use:
% [x_px,y_px] = transformPointsInverse(Tpx2rw,pos3D(:,1),pos3D(:,2),pos3D(:,3));
% ------------------------------------------------------------------------------------------

%% Definition of optional variables it they do not exist.
if ~exist('blackDots','var')
    blackDots=true;
end
if ~exist('extension','var')
    extension='tif';
end
if ~exist('FirstPlane','var')
    FirstPlane=1;
end
if ~exist('FirstCam','var')
    FirstCam=1;
end

% Total number of cameras
Ncam = numel(camName);
% Total number of calibration planes
NbzPlanes = numel(zPlanes);

% Initialisation
xyzRef(1).ref(1:NbzPlanes,:) = repmat([0 0 0],NbzPlanes,1); % pixels
xyzRef(2).ref(1:NbzPlanes,:) = repmat([0 0 0],NbzPlanes,1); % pixels
xyzRef(3).ref(1:NbzPlanes,:) = repmat([0 0 0],NbzPlanes,1); % pixels
xyzRef(4).ref(1:NbzPlanes,:) = repmat([0 0 0],NbzPlanes,1); % pixels

%% Let's treat every calibration pictures
for kz = FirstPlane:numel(zPlanes)
    z = zPlanes(kz)
    for kcam = FirstCam:Ncam
        filename = sprintf('%s/CalibrationPlan_%d_cam%d.%s',dirIn, kz, kcam, extension);        
        %% calib 2D detect white dots over dark background
        if blackDots
            Img=imcomplement(imread(filename));
        else
            Img=imread(filename);
        end
        ddisk = strel('disk',dotSize);
        ImBk = imopen(Img,ddisk);
        Img = imsubtract(Img,ImBk); % Background substraction
        
        % Treatment of the picture (by hand)
        [pimg,pos3D,T3px2rw,T3rw2px,T1rw2px,aRoi]=calib2D(Img,th,dotSize,gridSpace,lnoise,kz,z,xyzRef(kcam).ref(kz,:),kcam,dirIn);
    end
end

%% make calib structure
for kz = 1:numel(zPlanes)
    for kcam = 1:Ncam
        kz
        kcam
        load([dirIn filesep 'calib2D_' num2str(kz) '_cam' num2str(kcam) '.mat']);
        calib(kz,kcam).posPlane = zPlanes(kz);
        calib(kz,kcam).pimg = pimg;
        calib(kz,kcam).pos3D = pos3D;
        calib(kz,kcam).movedpoints = movedpoints;
        calib(kz,kcam).addedpoints = addedpoints;
%         calib(kz,kcam).T1rw2px = fitgeotrans(pimg,pos3D(:,1:2),'projective');
        calib(kz,kcam).T3rw2px = fitgeotrans(pimg,pos3D(:,1:2),'polynomial',3);
        calib(kz,kcam).T1px2rw = fitgeotrans(pos3D(:,1:2),pimg,'projective');
        calib(kz,kcam).T3px2rw = fitgeotrans(pos3D(:,1:2),pimg,'polynomial',3);
        calib(kz,kcam).cHull=convHullpimg;
        calib(kz,kcam).name = kcam;
        calib(kz,kcam).dirPlane=[1,2,3]; % Mire is displaced along axis 3, and 1,2 because the transformation provide x and y in real world

    end
end

save(sprintf('%s/calib.mat',dirIn),'calib');

end

