function BackgroundComputation(session,ManipName,CamNum,StartFrame,EndFrame ,Step,format)
%%% Compute the background for the pictures of the camera NumCam in the
%%% experiment ManipName
%%% Between picture StartFrame and EndFrame ( every Step frame) it takes
%%% the maximal/minimal/mean intensity for each pixel.
%----------------------------------------------------------------------------
%%% Parameters : 
%%%     session      : Path to the achitecture root (2 fields: session.input_path
% and session.output_path)
%%%     ManipName    : Name of the folder experiment
%%%     NumCam       : number of the camera studied
%%%     StartFrame   : number of first frame
%%%     EndFrame     : number of the last frame
%%%     Step         : step between 2 frame taken for computation
%%%     format (optional)     : picture names. By defaut it is '%05d.tif'.
%%%     The beginning of picture names has to be
%%%     %ManipName_cam%CamNum_%format.
%------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% By defaut format='%05d.tif'
if ~exist('format','var')
    format='%05d.tif';
end

%% Definitions of folders
folderin = sprintf("%sDATA/",session.input_path);
folderout = sprintf("%sProcessed_DATA/",session.output_path);

% If Step is not define take a step such as the average will be done over
% 1000 pictures.
if (~ exist("Step","var") && (EndFrame-StartFrame)>=1000)
    Step = floor((EndFrame-StartFrame)/1000);
elseif (~ exist("Step","var") && (EndFrame-StartFrame)<1000)
    Step = 1;
end

compteur = 0;
for kframe=StartFrame:Step:EndFrame
    ImgName = sprintf(join(['%s/%s_cam%d_' format],''),folderin, ManipName, CamNum, kframe);
    fprintf("%s \n",ImgName);
    
    if kframe==StartFrame
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
BackgroundMean = uint16(BackgroundMean/compteur);

save(sprintf("%s%s/Background_cam%d.mat",folderout,ManipName,CamNum),'BackgroundMax','BackgroundMean','BackgroundMin')


