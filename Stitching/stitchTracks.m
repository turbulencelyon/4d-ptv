
function StitchedTraj = stitchTracks(traj,dfmax,filepath,dxmax,dvmax,lmin)
% 2020 - David Dumont -> loop for df=1:dfmax
% 2017 - Mickaël Bourgoin
%
% Reconnect trajectories when their particles have been missed for <dfmax>
% frames maximum.
% ____________________________________________________________________________
% INPUTS
% traj   : structure composed of trajectory data
% dfmax  : maximum number of missing frames tolerated
% filepath
% dxmax  : maximum space difference tolerated between expected and real
% positions
% dvmax  : maximum relative velocity difference tolerated
% lmin   : minimum length to tru to reconnect a trajectory.
%
% OUTPUTS
% StitchedTraj(kff).fstitch         : index of interpolated frames
% StitchedTraj(kff).nbfstitch       : number of interpolated frames
% StitchedTraj(kff).dfstitch        : list of number of missing frames for each stitching
% StitchedTraj(kff).dXstitch        : list of dX value for each stitching
% StitchedTraj(kff).dVstitch        : list of dV value for each stitching
% StitchedTraj(kff).nbstitch        : number of stitch for kff trajectory
% StitchedTraj(kff).frames          : frame number
% StitchedTraj(kff).x               : x positions
% StitchedTraj(kff).y               : y positions
% StitchedTraj(kff).z               : z positions
% StitchedTraj(kff).L               : trajectory length
% StitchedTraj(kff).ntraj           : trajectory number in the tracking referential
% StitchedTraj(kff).nmatch          : index of each particle composing the trajectory. Index in the matching file.
% ____________________________________________________________________________

tic

%% We suppress all trajectories shorter than lmin to only treat longer ones
traj([traj.L]<lmin)=[];
Size = numel(traj);

%% new version (in construction) 7 x faster
L = vertcat(traj.L);
frames = vertcat(traj.frames);
x = vertcat(traj.x);
y = vertcat(traj.y);
z = vertcat(traj.z);

index_f = cumsum(L);
index_i = vertcat(1,index_f(1:end-1)+1);

fi=frames(index_i);
ff=frames(index_f);
xi=x(index_i);
xf=x(index_f);
yi=y(index_i);
yf=y(index_f);
zi=z(index_i);
zf=z(index_f);

index_vi = zeros(Size,floor(lmin/2)+1);
index_vf = zeros(Size,floor(lmin/2)+1);
for i=1:Size
   index_vi(i,:) = index_i(i):index_i(i)+floor(lmin/2);
   index_vf(i,:) = index_f(i)-floor(lmin/2):index_f(i);
end

vxi = mean(diff(x(index_vi),1,2),2);              % x velocity at the beginning of each trajectory (computed using lmin/2 first points)
vxf = mean(diff(x(index_vf),1,2),2);
vyi = mean(diff(y(index_vi),1,2),2);
vyf = mean(diff(y(index_vf),1,2),2);
vzi = mean(diff(z(index_vi),1,2),2);
vzf = mean(diff(z(index_vf),1,2),2);

%% Initialisation of stitching parameters
Nstitch = 0; 
StitchedTraj=traj;
if ~isfield(StitchedTraj,'nbfstitch')
    [StitchedTraj.nbfstitch] = deal(0);
end
if ~isfield(StitchedTraj,'fstitch')
    [StitchedTraj.fstitch] = deal([]);
end
if ~isfield(StitchedTraj,'dXstitch')
    [StitchedTraj.dXstitch] = deal([]);
end
if ~isfield(StitchedTraj,'dVstitch')
    [StitchedTraj.dVstitch] = deal([]);
end
if ~isfield(StitchedTraj,'dfstitch')
    [StitchedTraj.dfstitch] = deal([]);
end
if ~isfield(StitchedTraj,'nbstitch')
    [StitchedTraj.nbstitch] = deal(0);
end

%% Loop over ff (frame number at the end of each trajectory
fprintf("Let's start loop over trajectories\n")
kff = 1;
while (kff <= numel(ff))
    if rem(kff,1000)==0
        fprintf('Trajectories up to %d reconnected \n',kff)
    end
    %% Are there some trajectories which start in the frame where the kff trajectory ends?
    ii = find(and(fi>=ff(kff)+2,fi<=ff(kff)+dfmax+1)); % fi>=ff+2 and fi<=ff+dfmax     (ff | fi) at the minimum (| means missing frame) and (ff dfmax missing frames fi)
 
    % Yes
    if ~isempty(ii)
        df = fi(ii)-ff(kff);                                                                                                  % Trajectories distance in terms of frame number
        %dV = sqrt((vxi(ii)-vxf(kff)).^2+(vyi(ii)-vyf(kff)).^2+(vzi(ii)-vzf(kff)).^2)/sqrt(vxf(kff)^2+vyf(kff)^2+vzf(kff)^2);  % Relative velocity norm difference dV = sqrt(vi-vf)/sqrt(vf)
        dVx=(vxi(ii)-vxf(ii))/vxf(kff);
        dVy=(vyi(ii)-vyf(ii))/vyf(kff);
        dVz=(vzi(ii)-vzf(ii))/vzf(kff);
        dX = sqrt((xi(ii)-xf(kff)-df*vxf(kff)).^2+(yi(ii)-yf(kff)-df*vyf(kff)).^2+(zi(ii)-zf(kff)-df*vzf(kff)).^2);           % Position norm difference dX = sqrt(Xi-Xf-df*vf)
        
        %% Is there any trajectory close spatially with a similar velocity?
        Imatch=find((dX<dxmax) & (dVx<dvmax) & (dVy<dvmax) & (dVz<dvmax));
        % Yes
        if ~isempty(Imatch)
            if numel(Imatch)>1 % If there are several possible match we take the closest one.
                [~, imin] = min(dX(Imatch)+sqrt(dVx(Imatch).^2+dVy(Imatch).^2+dVz(Imatch).^2)); % Pas homogène !!!
                Imatch = Imatch(imin);
                imatch = ii(Imatch);
            else
                imatch = ii(Imatch); % Imatch is the index (in ii) among trajectories which start close in time to ff(kff) and imatch is the index in among all trajectories.
            end
            %% Construction of data in the range where particle is missing: we interpolate data.
            if df(Imatch) > 1
                finterp = ( ff(kff)+1 :ff(kff) + df(Imatch)-1 )';
                xinterp = interp1([StitchedTraj(kff).frames(end-floor(lmin/2):end); StitchedTraj(imatch).frames(1:floor(lmin/2)+1)],[StitchedTraj(kff).x(end-floor(lmin/2):end); StitchedTraj(imatch).x(1:floor(lmin/2)+1)],(ff(kff)+1 : ff(kff)+df(Imatch)-1) )';
                yinterp = interp1([StitchedTraj(kff).frames(end-floor(lmin/2):end); StitchedTraj(imatch).frames(1:floor(lmin/2)+1)],[StitchedTraj(kff).y(end-floor(lmin/2):end); StitchedTraj(imatch).y(1:floor(lmin/2)+1)],(ff(kff)+1 : ff(kff)+df(Imatch)-1) )';
                zinterp = interp1([StitchedTraj(kff).frames(end-floor(lmin/2):end); StitchedTraj(imatch).frames(1:floor(lmin/2)+1)],[StitchedTraj(kff).z(end-floor(lmin/2):end); StitchedTraj(imatch).z(1:floor(lmin/2)+1)],(ff(kff)+1 : ff(kff)+df(Imatch)-1) )';
            else  % a priori impossible to come there (see l.90) but let's be farsighted.
                finterp = [];
                xinterp = [];
                yinterp = [];
                zinterp = [];
            end

            StitchedTraj(kff).fstitch = [StitchedTraj(kff).fstitch; (ff(kff):fi(imatch))'; StitchedTraj(imatch).fstitch];         % The first and last frames here belong to previous and following trajectories. Conserved to get something when there are not missing frames
            StitchedTraj(kff).nbfstitch = numel(StitchedTraj(kff).fstitch);
            StitchedTraj(kff).frames = [StitchedTraj(kff).frames;  finterp; StitchedTraj(imatch).frames];
            StitchedTraj(kff).x = [StitchedTraj(kff).x;  xinterp; StitchedTraj(imatch).x];
            StitchedTraj(kff).y = [StitchedTraj(kff).y;  yinterp; StitchedTraj(imatch).y];
            StitchedTraj(kff).z = [StitchedTraj(kff).z;  zinterp; StitchedTraj(imatch).z];
            StitchedTraj(kff).L = StitchedTraj(kff).L + (fi(imatch)-ff(kff)-1) + StitchedTraj(imatch).L;
            StitchedTraj(kff).nmatch = [StitchedTraj(kff).nmatch; nan((fi(imatch)-ff(kff)-1),1); StitchedTraj(imatch).nmatch];   % We put some NaN where data are interpolated between two stitched trajectories
            StitchedTraj(kff).dfstitch = [StitchedTraj(kff).dfstitch; df(Imatch)];
            StitchedTraj(kff).dXstitch = [StitchedTraj(kff).dXstitch; dX(Imatch)];
            StitchedTraj(kff).dVstitch = [StitchedTraj(kff).dVstitch; [dVx(Imatch),dVy(Imatch),dVz(Imatch)]'];
            StitchedTraj(kff).nbstitch = numel(StitchedTraj(kff).dfstitch);
            
            % Actualisation of end characteristics of kff trajectory
            ff(kff) = ff(imatch);
            xf(kff) = xf(imatch);
            vxf(kff) = vxf(imatch);
            yf(kff) = yf(imatch);
            vyf(kff) = vyf(imatch);
            zf(kff) = zf(imatch);
            vzf(kff) = vzf(imatch);

            % Suppression of imatch trajectory properties
            StitchedTraj(imatch) = [];
            ff(imatch)=[];
            fi(imatch)=[];
            xf(imatch)=[];
            yf(imatch)=[];
            zf(imatch)=[];
            vxf(imatch)=[];
            vyf(imatch)=[];
            vzf(imatch)=[];
            xi(imatch)=[];
            yi(imatch)=[];
            zi(imatch)=[];
            vxi(imatch)=[];
            vyi(imatch)=[];
            vzi(imatch)=[];

            Nstitch = Nstitch + 1;
            % We do not have to do kff = kff +1 because we have reconnected
            % kff trajectory so we have to look at its new end. It is not
            % possible to have imatch<kff because kff represents time. So
            % no indexing problem.
            
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
end

fprintf('\n%d tracks were stitched\n',Nstitch)
%% Save in .h5
if exist('lminSupr','var')
   StitchedTraj([StitchedTraj.L]<lminSupr)=[]; 
end
OutputFileName = fullfile([filepath '.h5'])
fprintf("StitchedtracksSides saved as %s\n", OutputFileName)
fields=fieldnames(StitchedTraj);
for k=1:length(fields) %for each field
    fieldk=fields{k};
    if length(vertcat(StitchedTraj.(fieldk)))<1000
        chunksize = length(vertcat(StitchedTraj.(fieldk))); % chunksize has to be smaller or equal to fieldk size. This case occurs for processing of small frame packets.
    else
        chunksize = 1000;
    end
    h5create(OutputFileName,['/' fieldk],[length(vertcat(StitchedTraj.(fieldk))) 1],'ChunkSize',[chunksize 1],'Deflate',9)
    h5write(OutputFileName,['/' fieldk],vertcat(StitchedTraj.(fieldk)))
end
% Saving of stitching parameters
for param=["dfmax","dxmax","dvmax","lmin"]
    h5create(OutputFileName,sprintf('/stitchingparameters/%s',param),1)
    h5write(OutputFileName,sprintf('/stitchingparameters/%s',param),eval(param))
end

elapsed = toc;
fprintf("elapsed time = %d s",elapsed)

end
