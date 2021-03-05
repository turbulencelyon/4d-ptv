function CC = CenterFinding2D(session,ManipName,CamNum,firstFrame,nframes,th,sz,Test,BackgroundType,format)
%%% Detect particles position in picture and provides their positions in
%%% px.
%--------------------------------------------------------------------------------
%%% Parameters :
%%%     session                    : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
%%%     ManipName                  : Name of the folder experiment
%%%     NumCam                     : number of the camera studied
%%%     firstFrame                 : number of the first frame
%%%     nframes                    : total number of pictures
%%%     th                         : threshold
%%%     sz                         : typical size of the particles
%%%     Test                       : true-> test mode, false-> classic mode (optional)
%%%     BackgroundType (optional)  : determine which background is substracted to pictures. By defaut is equal to BackgroundMean,
%%%     format (optional)          : picture names. By defaut it is '%05d.tif'.
%%%     The beginning of picture names has to be %ManipName_cam%CamNum_%format
%--------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

%% Test if Test exist or not
if ~exist('Test','var')
    Test=false;
end

%% Test if BackgroundType exist or not
if ~exist('BackgroundType','var')
    BackgroundType="BackgroundMean";
end

% By defaut format='%05d.tif'
if ~exist('format','var')
    format='%05d.tif';
end

%% Definition of folders
fprintf(ManipName);
folderin = sprintf("%sDATA/%s/cam%d",session.input_path,ManipName,CamNum)
folderout = sprintf("%sProcessed_DATA/%s",session.output_path,ManipName)
BackgroundFile = sprintf("%s/Background_cam%d.mat",folderout,CamNum);

%% Find centers
if exist(folderout,'dir')==0
    mkdir(folderout);
end

load(BackgroundFile,'BackgroundMin','BackgroundMax','BackgroundMean')
%% Choice of background type
if BackgroundType=="BackgroundMean"
    Background=BackgroundMean;
elseif BackgroundType=="BackgroundMax"
    Background=BackgroundMax;
elseif BackgroundType=="BackgroundMin"
    Background=BackgroundMin;
end

if ~Test
    for kframe=firstFrame:nframes
        disp(kframe);
        ImgName = sprintf(join(['%s/%s_cam%d_' format],''),folderin, ManipName, CamNum, kframe);
        fprintf("%s \n",ImgName);
        Im = imread(ImgName) - Background;
        Nx = size(Im,2);
        Ny = size(Im,1);

        out=pkfnd(Im,th,sz); % Provides intensity maxima positions
        npar = size(out,1);
        
        %% We keep only spots with a gaussian shape
        cnt = 0;
        x = [];
        y = [];
        for j = 1:npar
            Nwidth = 1;
            if (out(j,2)-Nwidth >0)&&(out(j,1)-Nwidth>0)&&(out(j,2)+Nwidth<Ny)&&(out(j,1)+Nwidth<Nx)
                cnt = cnt+1;

                Ip = double(Im(out(j,2)-Nwidth:out(j,2)+Nwidth,out(j,1)-Nwidth:out(j,1)+Nwidth));

                x(end+1) = out(j, 1) + 0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
                y(end+1) = out(j, 2) + 0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));
            end
        end

        CC(kframe).X=x;
        CC(kframe).Y=y;
    end
    
    %% Centers saving into a .mat file
    save(sprintf("%s/centers_cam%d.mat",folderout,CamNum),"CC",'nframes','-v7.3')
else
    kframe=1
    ImgName = sprintf(join(['%s/%s_cam%d_' format],''),folderin, ManipName, CamNum, kframe);
    fprintf("%s \n",ImgName);
    Im = imread(ImgName) - Background;
    figure("NumberTitle","Off","Name",sprintf("RAW picture, cam %d, frame %d",CamNum,kframe))
    imshow(imread(ImgName))
    colormap gray
    figure("NumberTitle","Off","Name",sprintf("%s, cam %d",BackgroundType,CamNum))
    imshow(BackgroundMin)
    colormap gray
    figure("NumberTitle","Off","Name",sprintf("RAW picture - Background, cam %d, frame %d",CamNum,kframe))
    imshow(Im)
    colormap gray
    colorbar
    
%     if exist('Erosion','var')
%         se = strel('disk',1);
%         Imerode = imerode(Im,se);
%         Imdilate = imdilate(Imerode,se);
%         figure()
%         imagesc(Imdilate)
%         axis image
%         colormap gray
% 
%         Image = Im;
%         Im = Imdilate;
%     end

    %% Tracé de l'histogramme des intensités pour définir le seuil
    fig = figure('NumberTitle','Off','Name','Intensity histogram');
    histogram(Im,1000)
    xlabel("Intensity")
    ylabel("Number")
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    
    Nx = size(Im,2);
    Ny = size(Im,1);

    out=pkfnd(Im,th,sz); % Provides intensity maxima positions
    npar = size(out,1);
        
    %% We keep only spots with a gaussian shape
    cnt = 0;
    x = [];
    y = [];
    for j = 1:npar
        Nwidth = 1;
        if (out(j,2)-Nwidth >0)&&(out(j,1)-Nwidth>0)&&(out(j,2)+Nwidth<Ny)&&(out(j,1)+Nwidth<Nx)
            cnt = cnt+1;

            Ip = double(Im(out(j,2)-Nwidth:out(j,2)+Nwidth,out(j,1)-Nwidth:out(j,1)+Nwidth));

            x(end+1) = out(j, 1) + 0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
            y(end+1) = out(j, 2) + 0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));
        end
    end
    CC(kframe).X=x;
    CC(kframe).Y=y;

    fprintf("%d treated \n",kframe)

    %% Let's plot picture and detected points on a graph !!! Be careful the vertical axis is reversed compared to reality !!!
    figure('NumberTitle','Off','Name',sprintf("frame %d, %d detected points",kframe,numel(x)))
    imshow(Im)
    colormap gray
    
    hold on
    plot(flip(x),flip(y),'r+')
end

