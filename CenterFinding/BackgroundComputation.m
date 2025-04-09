function BackgroundComputation(session,ManipName,CamNum,firstFrame,endFrame,Step,format)
%%% Compute the background for the pictures of the camera NumCam in the
%%% experiment ManipName
%%% Between picture firstFrame and endFrame (every Step frame) it takes
%%% the maximal/minimal/mean intensity for each pixel.
%----------------------------------------------------------------------------
%%% Parameters : 
%%%     session      : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
%%%     ManipName    : Name of the folder experiment
%%%     NumCam       : number of the camera studied
%%%     firstFrame   : number of first frame
%%%     endFrame     : number of the last frame
%%%     Step         : step between 2 frame taken for computation
%%%     format (optional)     : picture names. By defaut it is '%05d.tif'.
%%%     The beginning of picture names has to be
%%%     %ManipName_cam%CamNum_%format.
%------------------------------------------------------------------------------
% 2020-2021 : D. Dumont (adapted from M. Bourgoin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% By defaut format='%05d.tif'
if ~exist('format','var')
    format='%05d.tif';
end

%% Definitions of folders
folderin = fullfile(session.input_path, 'DATA', ManipName, ['cam' num2str(CamNum)])
folderout = fullfile(session.output_path, 'Processed_DATA', ManipName)

%% Creation of output Folder if it does not exist
if ~isfolder(folderout) 
    mkdir(folderout)
end

% If Step is not define take a step such as the average will be done over
% 1000 pictures.
if (~ exist("Step","var") && (endFrame-firstFrame)>=1000)
    Step = floor((endFrame-firstFrame)/1000);
elseif (~ exist("Step","var") && (endFrame-firstFrame)<1000)
    Step = 1;
end

BaseName = join([ManipName '_cam' num2str(CamNum) '_' format],'');
compteur = 0;
for kframe=firstFrame:Step:endFrame
    ImgName = fullfile(folderin,sprintf(BaseName, kframe));
    fprintf("%s \n",ImgName);
    
    if kframe==firstFrame
        BackgroundMax = imread(ImgName);
        BackgroundMean = uint32(imread(ImgName));
        BackgroundMin = imread(ImgName);
    else
        BackgroundMax = max(cat(3,BackgroundMax,imread(ImgName)),[],3);
        BackgroundMean = BackgroundMean + uint32(imread(ImgName));
        BackgroundMin = min(cat(3,BackgroundMin,imread(ImgName)),[],3);
    end    
    compteur = compteur + 1;
end

%We save the background in the same format than the pictures
Picture = imread(ImgName);
if isa(Picture, 'uint8')
    BackgroundMean = uint8(BackgroundMean/compteur);
elseif isa(Picture, 'uint16')
    BackgroundMean = uint16(BackgroundMean/compteur);
else
    fprintf("Pictures have to be either 8 or 16-bits");
end

save(fullfile(folderout,['Background_cam' num2str(CamNum) '.mat']),'BackgroundMax','BackgroundMean','BackgroundMin')


