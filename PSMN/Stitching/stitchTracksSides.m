function StitchedTrajSides = stitchTracksSides(traj,minframe,maxframe,NbFramePerJobTracking,dfmax,dxmax,dvmax,lmin,filepath)
% 2020 - David Dumont -> loop for df=1:dfmax
% 2017 - Mickaël Bourgoin

% Reconnect trajectories comming from two differents files. Typically, if
% we have trajectories in frames 1 to 10 in <file1> and trajectories in
% frames 11 to 20 in <file2>, it stitches trajectories1 to trajectories2
% arround frame 10+/-dfmax exclusively. If dfmax=0, it is just like predictive
% tracking.
% ____________________________________________________________________________
% INPUTS
% traj                   : structure composed of trajectory data
% minframe               : number of the first frame to treat
% maxframe               : number of the last frame to treat
% NbFramePerJobTracking  : number of frame in tracking jobs
% dfmax                  : maximum number of missing frames tolerated 
% dxmax                  : maximum space difference tolerated between expected and real positions
% dvmax                  : maximum relative velocity difference tolerated
% lmin                   : minimum length to tru to reconnect a trajectory
% filepath (optional)    : if exist, save the output as filepath. If not exist, does not save output.
%
% OUTPUTS
% StitchedTrajSides(kff).fstitch         : index of interpolated frames
% StitchedTrajSides(kff).nbfstitch       : number of interpolated frames
% StitchedTrajSides(kff).dfstitch        : list of number of missing frames for each stitching
% StitchedTrajSides(kff).dXstitch        : list of dX value for each stitching
% StitchedTrajSides(kff).dVstitch        : list of dV value for each stitching
% StitchedTrajSides(kff).nbstitch        : number of stitch for kff trajectory
% StitchedTrajSides(kff).frames          : frame number
% StitchedTrajSides(kff).x               : x positions
% StitchedTrajSides(kff).y               : y positions
% StitchedTrajSides(kff).z               : z positions
% StitchedTrajSides(kff).L               : trajectory length
% StitchedTrajSides(kff).ntraj           : trajectory number in the tracking referential
% StitchedTrajSides(kff).nmatch          : index of each particle composing the trajectory. Index in the matching file.
% ____________________________________________________________________________

tic
%% We only treat  trajectories longer than lmin
traj([traj.L]<lmin)=[];
Size = numel(traj);

%% Let's select only trajectories which start/end at the frames at the end/start of file +/- dfmax
fi=zeros(Size,1);
ff=zeros(Size,1);
for i=1:Size
    fi(i) = traj(i).frames(1);        % frame number at the beginning of each trajectory
    ff(i) = traj(i).frames(end);      % frame number at the end of each trajectory
end

%% Reasearch of all trajectories ending between frames nframe-dfmax-1 and nframe-1 and trajectories starting between frames nframe and nframe+dfmax
for nframe=ListFileFrame[2:end-1] %(minframe+NbFramePerJobTracking):NbFramePerJobTracking:(maxframe-NbFramePerJobTracking+1)
    if nframe==minframe+NbFramePerJobTracking
        FI = and(fi>=nframe,fi<=nframe+dfmax);
        FF = and(ff>=nframe-dfmax-1,ff<nframe);
    else
        FI = or(FI, and(fi>=nframe,fi<=nframe+dfmax));
        FF = or(FF, and(ff>=nframe-dfmax-1,ff<nframe));
    end
end

%% Creation of trajectory indexes for the two types of trajectories
I_fi = find(FI);
I_ff = find(FF);
fprintf("%d trajectories starting close to file limits\n",numel(I_fi))
fprintf("%d trajectories ending close to file limits\n",numel(I_ff))

%% Creation of list containning frame number, position and velocity at the start/end of each starting/ending trajectory (I_fi/I_ff)
fprintf("Construction of pre-data...\n")
fi=zeros(numel(I_fi),1);
xi=zeros(numel(I_fi),1);
yi=zeros(numel(I_fi),1);
zi=zeros(numel(I_fi),1);
vxi=zeros(numel(I_fi),1);
vyi=zeros(numel(I_fi),1);
vzi=zeros(numel(I_fi),1);
% for starting trajectories
for i=1:numel(I_fi)
    fi(i) = traj(I_fi(i)).frames(1);                                    % frame number at the beginning of each trajectory
    xi(i) = traj(I_fi(i)).x(1);                                         % x position at the beginning of each trajectory
    vxi(i) = mean(diff(traj(I_fi(i)).x(1:floor(lmin/2)+1)));              % x velocity at the beginning of each trajectory (computed using lmin/2 first points)
    yi(i) = traj(I_fi(i)).y(1);                                         % y position at the beginning of each trajectory
    vyi(i) = mean(diff(traj(I_fi(i)).y(1:floor(lmin/2)+1)));              % y velocity at the beginning of each trajectory (computed using lmin/2 first points)
    zi(i) = traj(I_fi(i)).z(1);                                         % z position at the beginning of each trajectory
    vzi(i) = mean(diff(traj(I_fi(i)).z(1:floor(lmin/2)+1)));              % z velocity at the beginning of each trajectory (computed using lmin/2 first points)
end
% for ending trajectories
ff=zeros(numel(I_ff),1);
xf=zeros(numel(I_ff),1);
yf=zeros(numel(I_ff),1);
zf=zeros(numel(I_ff),1);
vxf=zeros(numel(I_ff),1);
vyf=zeros(numel(I_ff),1);
vzf=zeros(numel(I_ff),1);
for i=1:numel(I_ff)
    ff(i) = traj(I_ff(i)).frames(end);                                  % frame number at the end of each trajectory
    xf(i) = traj(I_ff(i)).x(end);                                       % x position at the end of each trajectory
    vxf(i) = mean(diff(traj(I_ff(i)).x(end-floor(lmin/2):end)));      % x velocity at the end of each trajectory (computed using lmin/2 last points)
    yf(i) = traj(I_ff(i)).y(end);                                       % y position at the end of each trajectory
    vyf(i) = mean(diff(traj(I_ff(i)).y(end-floor(lmin/2):end)));      % y velocity at the end of each trajectory (computed using lmin/2 last points)
    zf(i) = traj(I_ff(i)).z(end);                                       % z position at the end of each trajectory
    vzf(i) = mean(diff(traj(I_ff(i)).z(end-floor(lmin/2):end)));      % z velocity at the end of each trajectory (computed using lmin/2 last points)
end

%% Initialisation of stitching parameters
Nstitch = 0; 
StitchedTrajSides = traj;
if ~isfield(StitchedTrajSides,'nbfstitch')
    [StitchedTrajSides.nbfstitch] = deal(0);
end
if ~isfield(StitchedTrajSides,'fstitch')
    [StitchedTrajSides.fstitch] = deal([]);
end
if ~isfield(StitchedTrajSides,'dXstitch')
    [StitchedTrajSides.dXstitch] = deal([]);
end
if ~isfield(StitchedTrajSides,'dVstitch')
    [StitchedTrajSides.dVstitch] = deal([]);
end
if ~isfield(StitchedTrajSides,'dfstitch')
    [StitchedTrajSides.dfstitch] = deal([]);
end
if ~isfield(StitchedTrajSides,'nbstitch')
    [StitchedTrajSides.nbstitch] = deal(0);
end
[StitchedTrajSides.stitched] = deal(0);
%% Loop over ff (frame number at the end of each ending trajectory
fprintf("Let's start loop over trajectories\n")
kff = 1;
while (kff <= numel(ff))
    if rem(I_ff(kff),1000)==0
        fprintf('Trajectories up to %d reconnected \n',kff)
    end
    % Check if this trajectory was not already stitched
    if StitchedTrajSides(I_ff(kff)).stitched==0
        %% Are there some trajectories which start in the frame where the kff trajectory ends? Pay Attention it is different from stitchTracks because here fi>=ff+1 and not fi>=ff+2
        ii = find(and(fi>=ff(kff)+1,fi<=ff(kff)+dfmax+1)); % fi>=ff+1 and fi<=ff+dfmax+1     (ff fi) at the minimum (no missing frame) and (ff dfmax missing frames fi) at the maximum
        % Yes
        if ~isempty(ii)
            df = fi(ii)-ff(kff); % Trajectories distance in terms of frame number
            dVx=(vxi(ii)-vxf(kff))/vxf(kff);
            dVy=(vyi(ii)-vyf(kff))/vyf(kff);
            dVz=(vzi(ii)-vzf(kff))/vzf(kff);
            dX = sqrt((xi(ii)-xf(kff)-df*vxf(kff)).^2+(yi(ii)-yf(kff)-df*vyf(kff)).^2+(zi(ii)-zf(kff)-df*vzf(kff)).^2);           % Position norm difference dX = sqrt(Xi-Xf-df*vf)

            %% Is there any trajectory close spatially with a similar velocity?
            Imatch=find((dX<dxmax) & (abs(dVx)<dvmax) & (abs(dVy)<dvmax) & (abs(dVz)<dvmax));
            % Yes
            if ~isempty(Imatch)
                if numel(Imatch)>1 % If there are several possible match we take the closest one.
                    [~, imin] = min(dX(Imatch)+sqrt(dVx(Imatch).^2+dVy(Imatch).^2+dVz(Imatch).^2)); % Pas homogène !!!
                    Imatch = Imatch(imin);                      % Imatch is the index (in ii) among trajectories which start close in time to ff(kff)
                    imatch_fi = ii(Imatch);                     % imatch_fi is the index among I_fi
                    imatch = I_fi(ii(Imatch));                  % imatch if the index among all trajectories
                else
                    imatch_fi = ii(Imatch);                     % imatch_fi is the index among I_fi
                    imatch = I_fi(ii(Imatch));                  % imatch is the index among all trajectories
                end
                %% Construction of data in the range where particle is missing: we interpolate data.
                if df(Imatch) > 1
                    finterp = ( ff(kff)+1 :ff(kff) + df(Imatch)-1 )';
                    xinterp = interp1([StitchedTrajSides(I_ff(kff)).frames(end-floor(lmin/2):end); StitchedTrajSides(imatch).frames(1:floor(lmin/2)+1)],[StitchedTrajSides(I_ff(kff)).x(end-floor(lmin/2):end); StitchedTrajSides(imatch).x(1:floor(lmin/2)+1)],(ff(kff)+1 : ff(kff)+df(Imatch)-1) )';
                    yinterp = interp1([StitchedTrajSides(I_ff(kff)).frames(end-floor(lmin/2):end); StitchedTrajSides(imatch).frames(1:floor(lmin/2)+1)],[StitchedTrajSides(I_ff(kff)).y(end-floor(lmin/2):end); StitchedTrajSides(imatch).y(1:floor(lmin/2)+1)],(ff(kff)+1 : ff(kff)+df(Imatch)-1) )';
                    zinterp = interp1([StitchedTrajSides(I_ff(kff)).frames(end-floor(lmin/2):end); StitchedTrajSides(imatch).frames(1:floor(lmin/2)+1)],[StitchedTrajSides(I_ff(kff)).z(end-floor(lmin/2):end); StitchedTrajSides(imatch).z(1:floor(lmin/2)+1)],(ff(kff)+1 : ff(kff)+df(Imatch)-1) )';

                    else  % a priori impossible to come there (see l.90) but let's be farsighted.
                    finterp = [];
                    xinterp = [];
                    yinterp = [];
                    zinterp = [];
                end

                StitchedTrajSides(I_ff(kff)).fstitch = [StitchedTrajSides(I_ff(kff)).fstitch; (ff(kff):fi(imatch_fi))'; StitchedTrajSides(imatch).fstitch];                       % The first and last frames here belong to previous and following trajectories. Conserved to get something when there are not missing frames
                StitchedTrajSides(I_ff(kff)).nbfstitch = numel(StitchedTrajSides(I_ff(kff)).fstitch);
                StitchedTrajSides(I_ff(kff)).frames = [StitchedTrajSides(I_ff(kff)).frames;  finterp; StitchedTrajSides(imatch).frames];
                StitchedTrajSides(I_ff(kff)).x = [StitchedTrajSides(I_ff(kff)).x;  xinterp; StitchedTrajSides(imatch).x];
                StitchedTrajSides(I_ff(kff)).y = [StitchedTrajSides(I_ff(kff)).y;  yinterp; StitchedTrajSides(imatch).y];
                StitchedTrajSides(I_ff(kff)).z = [StitchedTrajSides(I_ff(kff)).z;  zinterp; StitchedTrajSides(imatch).z];
                StitchedTrajSides(I_ff(kff)).L = StitchedTrajSides(I_ff(kff)).L + (fi(imatch_fi)-ff(kff)-1) + StitchedTrajSides(imatch).L;
                StitchedTrajSides(I_ff(kff)).nmatch = [StitchedTrajSides(I_ff(kff)).nmatch; nan((fi(imatch_fi)-ff(kff)-1),1); StitchedTrajSides(imatch).nmatch];   % We put some NaN where data are interpolated between two stitched trajectories
                StitchedTrajSides(I_ff(kff)).dfstitch = [StitchedTrajSides(I_ff(kff)).dfstitch; df(Imatch)];
                StitchedTrajSides(I_ff(kff)).dXstitch = [StitchedTrajSides(I_ff(kff)).dXstitch; dX(Imatch)];
                StitchedTrajSides(kff).dVstitch = [StitchedTrajSides(kff).dVstitch; [dVx(Imatch),dVy(Imatch),dVz(Imatch)]'];
                StitchedTrajSides(I_ff(kff)).nbstitch = numel(StitchedTrajSides(I_ff(kff)).dfstitch);
                
                % Actualisation of end characteristics of kff trajectory
                ff(kff) = StitchedTrajSides(imatch).frames(end);
                xf(kff) = StitchedTrajSides(imatch).x(end);
                vxf(kff) = mean(diff(StitchedTrajSides(imatch).x(end-floor(lmin/2):end)));
                yf(kff) = StitchedTrajSides(imatch).y(end);
                vyf(kff) = mean(diff(StitchedTrajSides(imatch).y(end-floor(lmin/2):end)));
                zf(kff) = StitchedTrajSides(imatch).z(end);
                vzf(kff) = mean(diff(StitchedTrajSides(imatch).z(end-floor(lmin/2):end)));

                % Suppression of imatch trajectory properties
                StitchedTrajSides(imatch).stitched = 1;
%                 StitchedTrajSides(imatch).frames = nan;
                StitchedTrajSides(imatch).x = nan; % we kept a nan in x to remove it from dX or dV
%                 StitchedTrajSides(imatch).y = nan;
%                 StitchedTrajSides(imatch).z = nan;
%                 StitchedTrajSides(imatch).L = nan;
%                 StitchedTrajSides(imatch).nmatch = nan;
    %             ff(imatch_fi)=nan;
                fi(imatch_fi)=nan;
    %             xf(imatch_fi)=nan;
    %             yf(imatch_fi)=nan;
    %             zf(imatch_fi)=nan;
    %             vxf(imatch_fi)=nan;
    %             vyf(imatch_fi)=nan;
    %             vzf(imatch_fi)=nan;
                xi(imatch_fi)=nan;
                yi(imatch_fi)=nan;
                zi(imatch_fi)=nan;
                vxi(imatch_fi)=nan;
                vyi(imatch_fi)=nan;
                vzi(imatch_fi)=nan;
                I_fi(imatch_fi)=nan;

                Nstitch = Nstitch + 1;

            % No trajectory enough close spatially and with a similar velocity
            % so we jump to the following trajectory.
            else
                kff = kff + 1;
            end
        % No trajectories which start in the vicinity of ff(kff). So we jump to
        % the following trajectory.
        else
            kff = kff + 1;
        end
    % Trajectory already stitched    
    else
        kff = kff + 1;
    end
end

% Suppression of stitched trajectories
StitchedTrajSides([StitchedTrajSides.stitched]==1)=[];
StitchedTrajSides=rmfield(StitchedTrajSides,'stitched'); % supression of stitched field which is no more useful

fprintf('\n%d tracks were stitched\n',Nstitch)

%% Save in .h5
if exist('filepath','var')
    OutputFileName = sprintf('%s.h5',filepath);
    fields=fieldnames(StitchedTrajSides);
    for k=1:length(fields) %for each field
        fieldk=fields{k};
        if length(vertcat(StitchedTrajSides.(fieldk)))<1000
            chunksize = length(vertcat(StitchedTrajSides.(fieldk))); % chunksize has to be smaller or equal to fieldk size. This case occurs for processing of small frame packets.
        else
            chunksize = 1000;
        end
        h5create(OutputFileName,['/' fieldk],[length(vertcat(StitchedTrajSides.(fieldk))) 1],'ChunkSize',[chunksize 1],'Deflate',9)
        h5write(OutputFileName,['/' fieldk],vertcat(StitchedTrajSides.(fieldk)))
    end
    fprintf("StitchedtracksSides saved as %s\n", OutputFileName)
end
elapsed = toc;
fprintf("elapsed time = %d s\n",elapsed)

end