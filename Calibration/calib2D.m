function [pimg,pos3D,T3px2rw,T3rw2px,T1rw2px,aRoi]=calib2Df(Iimg,threshold,Size_dot,...
    gridspace,lnoise,kz,zPlane,xyzRef,camNumber,savepath,aRoi)
%
%[pimg,pos3D,T,Tinv,aRoi]=calib2Df(Iimg,threshold,Size_dot,gridspace,lnoise,zPlane,xyzRef,camNumber,savepath,aRoi)
%
% Takes an image of a calibration mask, applies some filtering in order to
% easily find centers of points on mask. Then asks for selection of a
% region of interest (IMPORTANT: start with upper-left corner and continue
% clockwise to avoid rotation of image). Then detects points inside this
% ROI, using the ROI contours to project the region into a rectangle, and
% thus sort the array of points. Once the points are sorted, knowing the
% real separation between points on grid, and the position of the plane in
% the z-coordinate, it gives the forward and inverse geometric
% transformations to the real-world coordinates on the plane, Tpx2rw (T)
% and Trw2px (Tinv) respectively.
%
%__________________________________________________________________________
%
% INPUTS
% Iimg      : original image
% threshold :
% Size_dot  : typical diameter, in pixels, of dots in mask image
% gridspace : grid spacing of calibration mask in mm
% lnoise    : typical lengthscale of image noise
% kz        : plan number (for file name)
% zPlane    : z-position of calibration mask plane, perpendicular to cam
% camNumber : number of camera (1, 2, 3, etc.. just for file saving)
% savepath  : path to directory where results will be saved
% aRoi      : [x_NorthWest y_NW ; x_NE y_NE ; x_SE y_SE ; x_SW y_SW]
%
% OUTPUTS
% pimg      : center coordinates in original image [in pixels]
% pos3D     : center coordinates in real world [in mm]
% T3rw2px   : transformation from real world to image cubic
% T1rw2px   : transformation from real world to image linear
% T3px2rw   : transformation from image to real world cubic
% aRoi      : coordinates of region of interest, in case we don't want to
%             re-run script without selecting it again manually.
%
% IMPORTANT! (the following is kind of confusing)
%
% To transform image to real world use:
% [x_rw,y_rw,z_rw] = transformPointsInverse(Trw2px,posimg(:,1),posimg(:,2));
%
% To transform real world positions to image positions use:
% [x_px,y_px] = transformPointsInverse(Tpx2rw,pos3D(:,1),pos3D(:,2),pos3D(:,3));
%
%__________________________________________________________________________

fname = ['plane_' num2str(zPlane) '_C' num2str(camNumber,'%03d')];

if ~exist('lnoise','var')
    lnoise = 1;
end

ncams         = 1;           % nmbr of cams in this calibration
movingaxis    = 'z';         % direction of displacement of plane
maskaxis      = ['x', 'y'];  % horizontal and vertical axes of plane
img_inversion = 'n';         % invert image colors [y/n]
roiConvHull   = 1 ;           % if 1 calibrate only within the ROI (Region Of Interest)
axispos = zPlane;            % position of plane in moving axis coordinate
pimg    = [];
pos3D   = [];

% invert image colors
if lower(img_inversion(1)) == 'y' 
    Iimg=imcomplement(Iimg);
end

% bandpass filters image according to level of noise and size of dots
kg = fspecial('gaussian',2*Size_dot,lnoise);
Iimg=imfilter(Iimg,kg);

[Npix_y, Npix_x] = size(Iimg);
% f = figure(nfig);
% p = uipanel('Parent',f,'BorderType','none'); 
% p.Title = ['zPlane = ' num2str(zPlane) '  ;  cam# = ' num2str(camNumber) '  ;  xyzRef = [' num2str(xyzRef) ']' ];
% p.TitlePosition = 'centertop'; 
% p.FontSize = 36;
% p.FontWeight = 'bold';
% sgtitle(['zPlane = ' num2str(zPlane) '  ;  cam# = ' num2str(camNumber) '  ;  xyzRef = [' num2str(xyzRef) ']']);
redo = 'y';

while strcmpi(redo, 'y')
    
    xc = []; yc = []; Ap = [];
    h1 = subplot('position', [0.02 0.3 0.37 0.6]);
    h2 = subplot('position', [0.02 0.08 0.37 0.15]);
    h3 = subplot('position', [0.4 0.08 0.59 0.8]);
    subplot(h1);
    imagesc(Iimg); %colormap(gray);
    set(gcf, 'Position', get(0, 'Screensize')); hold on;
    title(fname);
    
    %%
    %disp('Now please indicate a square or rectangular region of interest in real space');
    title('Now please indicate a square or rectangular region of interest in real space');
    
    if ~exist('aRoi','var')  % select Region of Interest for calibration
        [BW_roi,xi,yi] = roipoly;
        aRoi=[xi(1:4) yi(1:4)];
    else                     % load Region of Interest
        [BW_roi,xi,yi]=roipoly(Iimg,aRoi(1:4,1),aRoi(1:4,2));
        aRoi=[xi(1:4) yi(1:4)];
    end
    
    % use preliminary transform to project ROI onto a square and reshape
    % array of points into a squared array, in this way, its easier to
    % sort the centroid positions
    a0 = [0 0; 1 0; 1 1; 0 1];
    T0 = cp2tform(aRoi,a0,'projective');
%     T0 = fitgeotrans(aRoi,a0,'projective');
    
    %% Histogram Intensity plot
    Nhist = hist(reshape(double(Iimg), size(Iimg,1)*size(Iimg,2), 1), [0:255]);
    subplot(h2);
    semilogy(0:255, Nhist, 'b-');
    axis([0 255 1 10000]);
    th = threshold;
    
    % suppress structures connected to image border
    Ith = imclearborder(Iimg>th);
    
    % find structures connected components in binary image
    if roiConvHull == 1
        SS    = regionprops(logical(Ith.*BW_roi),'Centroid','Area');
    else
        SS    = regionprops(logical(Ith),'Centroid','Area');
    end
    Apsub = [SS.Area]';
    
    % discard small dots
    II        = find(Apsub<0.5*pi*Size_dot^2/4);
    SS(II)    = [];
    Apsub(II) = [];
    XX        = reshape([SS.Centroid],2,[])';
    xc        = XX(:,1);                      % x-position of centroids
    yc        = XX(:,2);                      % y-position of centroids
    
    if roiConvHull ==1
        % remove all points outside of ROI
        i_OutROI     = ~inpolygon(xc,yc,xi,yi);
        xc(i_OutROI) = [];
        yc(i_OutROI) = [];
    end
    
    Ap = [Ap ; Apsub];
    
    %% Image + calibration point plot
    subplot(h3);

    % transform image using projective transform generated with ROI
    [Itht,xx,yy] = imtransform(Iimg,T0,'Size',size(Ith));

    if img_inversion==0
        imagesc(xx,yy,Itht); colormap('gray');
%         imagesc(Itht); colormap('gray');
    else
        imagesc(xx,yy,imcomplement(Itht)); colormap('gray');
    end
    hold on;

    % transform centroids positions using the same projective transform
    xc0     = xc;
    yc0     = yc;
    [xc,yc] = tformfwd(T0,xc0,yc0);
    
    % plot centroids
    plot(xc, yc, 'r+');
    hold off
    redo = 'n';
end

Np = length(xc); % number of detected dots

% % manually remove/add/displace some points if it's more convenient
ind = ones(Np,1);
nrm = 0;
title({'Please click the particle centers that you want to remove. Right click the mouse when you are done.','defaut: suppress mode, a : adding mode, d: displacement mode'});
but = 1;
addedpoints = zeros(Np,1); % list of added points
movedpoints = zeros(Np,1); % list of moved points
hold on
while but ~= 3  % si on fait un clic droit alors but = 3
    [xrm, yrm, but] = ginput(1);
    if but == 1
        dist = (xc-xrm).^2+(yc-yrm).^2;
        [mindist, irm] = min(dist);
        % If the point is not yet selected, then it is removed
        if ind(irm) == 1
            plot(xc(irm), yc(irm), 'ro');
            ind(irm) = 0;
            nrm = nrm+1;
        else % if the point is already selected, then it is unselected 
            plot(xc(irm), yc(irm), 'wo'); % to remove red circle
            plot(xc(irm), yc(irm), 'r');
            ind(irm) = 1;
            nrm = nrm-1;
        end
    end
    % switch between supress (s) (defaut)/deplacement (d)/add (a) mode
    keyboard=get(gcf,'CurrentCharacter');
    if keyboard~='0' % has it changed from the dummy character?
    set(gcf,'CurrentCharacter','0'); % reset the character
        % deplacement mode : to deplace points
        if keyboard=='d'
            [xd, yd, but] = ginput(1);
            if but == 1
                dist = (xc-xd).^2+(yc-yd).^2;
                [mindist, id] = min(dist);
                plot(xc(id), yc(id), 'go'); % point becomes green to show which one is selected
                [xdn, ydn, but] = ginput(1);
                xc(id)=xdn;
                yc(id)=ydn;
                plot(xc(id), yc(id), 'ko'); % new point appears black
                plot(xc(id), yc(id), 'r+');
                movedpoints(id) = 1; % update of moved points list
            end
        end
        % adding mode : to add missing points
        if keyboard=='a'
           [xa, ya, but] = ginput(1); 
           xc(end+1)=xa;
           yc(end+1)=ya;
           plot(xc(end),yc(end),'co') % new point appears cyan
           plot(xc(end),yc(end),'r+')
           addedpoints(end+1) = 1; % update of added points list
           ind(end+1) = 1; % update of logical points list (which determine later if the point is conserved (1) or not (0) )
        end
    end
end
hold off

% expansion of movedpoints to get the same size than addedpoints
Np = length(xc); % size of xc before removing suppressed points
% If some points were added we expanse xc0 and yc0 to take it into account
if sum(addedpoints)>0
    xc0(Np) = 0; % the Np item is set to 0
    yc0(Np) = 0;
    SizeMovedPoints = size(movedpoints);
    if SizeMovedPoints(1)<Np % That will be true if the last added point was not moved. If it was moved, movedpoints would have already Np size
        movedpoints(Np)=0;
    end
end

%% Added and displaced points are known in the real world and in the rescaled image but not yet in the original image coordinates. 
% We applied an inverse transformation to get it (xc0, yc0 from xc,yc)
[xcad, ycad] = tforminv(T0, xc(logical(addedpoints+movedpoints)), yc(logical(addedpoints+movedpoints))); % only moved/added points are transformed

xc0(logical(addedpoints+movedpoints))=xcad;
yc0(logical(addedpoints+movedpoints))=ycad;

xc = xc(logical(ind));
yc = yc(logical(ind));
xc0 = xc0(logical(ind));
yc0 = yc0(logical(ind));

% Ap = Ap(logical(ind),:);
Np = length(xc);
disp(strcat(num2str(nrm), ' points have been removed. Please check the image again'));
fprintf('%d points have been moved \n',sum(movedpoints))
fprintf('%d points have been added \n', sum(addedpoints))

if img_inversion==0
    imagesc(xx,yy,Itht);
else
    imagesc(xx,yy,imcomplement(Itht));
end
hold on;


plot(xc, yc, 'r+');
hold off

% select base vectors in correct order (center, right and then top)
title('Now please indicate the three base point on the mask by click mouse on the thresholded image in order [0 0], [1 0], [0 1].');

% select and identify first base point: center
but = 0;
while but ~= 1
    [x0 y0 but] = ginput(1);
end
dist = (xc-x0).^2+(yc-y0).^2;
[mindist i0] = min(dist);
subplot(h3);
hold on
plot(xc(i0), yc(i0), 'bo');

i0xind = 0;
i0yind = 0;

% select and identify second base point: east 
but = 0;
while but ~= 1
    [x1 y1 but] = ginput(1);
end
dist = (xc-x1).^2+(yc-y1).^2;
[mindist i1] = min(dist);
plot(xc(i1), yc(i1), 'bo');

i1xind = 1;
i1yind = 0;

% select and identify third base point: north
but = 0;
while but ~= 1
    [x2 y2 but] = ginput(1);
end
dist = (xc-x2).^2+(yc-y2).^2;
[mindist i2] = min(dist);
plot(xc(i2), yc(i2), 'bo');
i2xind = 0;
i2yind = 1;

% Now determine the point coordinates
% first, form two base vectors on the mask
e1 = [i1xind-i0xind, i1yind-i0yind];
e2 = [i2xind-i0xind, i2yind-i0yind];

% The projection of these two vectors on image plane
e1p = [xc(i1)-xc(i0), yc(i1)-yc(i0)];
e2p = [xc(i2)-xc(i0), yc(i2)-yc(i0)];

e1pnorm = sum(e1p.^2);   % their squared norms
e2pnorm = sum(e2p.^2);

e1pe2p = sum(e1p.*e2p);  % their dot product

d = (e1pnorm*e2pnorm - e1pe2p*e1pe2p);		% the denominator

% calculate the coords of all points using the two base vectors
pind = zeros(Np, 2);
for i = 1:Np
    c = [xc(i)-xc(i0), yc(i)-yc(i0)]; % coord in pixel units
    A = (sum(c.*e1p)*e2pnorm - sum(c.*e2p)*e1pe2p)/d;
    B = (sum(c.*e2p)*e1pnorm - sum(c.*e1p)*e1pe2p)/d;
    pind(i,:) = [A B];
end

% Now calculate the two components of dots' 3D coordinates on the mask plane
pmask = zeros(Np,2);
originalpmask = zeros(Np,2);

for i = 1:Np
    originalpmask(i,1:2) = (e1*pind(i,1))+(e2*pind(i,2)) + [i0xind i0yind];
    pmask(i,1:2)         = round(originalpmask(i,1:2));
    NN(i)                = norm(pmask(i,1:2)-originalpmask(i,1:2));

%% for single plane target
     pmask(i,3)=0;
end

Np = size(pmask,1);

% check to see if there is any inconsistency
ncoll = 0;
icoll = [];

for i = 1:Np
    str = strcat('(',num2str(pmask(i,1)),',',num2str(pmask(i,2)),',',num2str(pmask(i,3)),')');
	plot(xc(i),yc(i),'yx');
    text(xc(i), yc(i), str, 'Color', 'm');
end

% check if there are repeated mask points
for i = 2:Np
    for j = 1:i-1
        if (pmask(i,1) == pmask(j,1)) && (pmask(i,2) == pmask(j,2))
            ncoll = ncoll+1;
            icoll = [icoll; [i j]];
        end
    end
end

if ncoll == 0
    disp('No conflicts found, but please still check particle coordinates.');
    title('No conflicts found, but please still check particle coordinates.');
else
    str = sprintf('Some particle coordinates are probably wrong. %d conflicts found: \n', ncoll);
    for  ic = 1:ncoll
        str = strcat(str, sprintf('(%d, %d)\n', pmask(icoll(ic),:)));
    end
    disp(str);
    title(str);
end

% We avait for Enter to valid the calibration
pause


pmask      = pmask + repmat(xyzRef,length(pmask),1);
pos3D      = pmask*gridspace;
pos3D(:,3) = pos3D(:,3)+ zPlane;   % position of dots in 3d [mm]
pimg       = [xc0 yc0];            % position of dots in image [pix]
pos2D      = pos3D(:,1:2);         % position of dots in 2d [mm]
movedpoints= movedpoints(logical(ind));
addedpoints= addedpoints(logical(ind));

% compute 3rd order polynomial spatial transformation from image points
% [in pixels] to 2d position in real-space on plane [in mm]
ttype = 'polynomial';
T3rw2px  = fitgeotrans(pimg,pos2D,ttype,3); % inverse transform
T3px2rw  = fitgeotrans(pos2D,pimg,ttype,3); % forward transform

ttype = 'projective';
T1rw2px  = fitgeotrans(pimg,pos2D,ttype); % inverse transform

% save results to file
convHullpimg = convhull(pimg(:,1),pimg(:,2));
save(sprintf('%s/calib2D_%d_cam%d',savepath, kz ,camNumber),...
    'zPlane','T3rw2px','T3px2rw','T1rw2px','pos3D','pimg','convHullpimg','movedpoints','addedpoints');
close
