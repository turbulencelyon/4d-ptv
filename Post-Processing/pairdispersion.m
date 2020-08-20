function delta = pairdispersion(traj)
% 05/2020 - David Dumont (inspirated by Thomas Basset's code
% Compute for all timestamp distances between each couple of trajectories.
% -------------------------------------------------------------------------------
% INPUT
%     traj : structure containing all trajectories with at least the following fields:
%     traj.x       x position
%     traj.y
%     traj.z 
%     traj.frames  frames number
%     traj.ntraj   trajectory index
%     traj.L       trajectory length
%     
% OUTPUT
%     delta : structure containing for each couple of trajectories
%     delta.dx
%     delta.dy
%     delta.dz   z distance...
%     delta.d    distance between the two trajectories
%     delta.dt   frame number (time)
%     delta.dn   trajectories indexes
%     delta.L    length
% -------------------------------------------------------------------------------

tic

% Let's start by supressing small trajectories (smaller than 75 frames)
traj([traj.L]<75)=[];

frames=vertcat(traj.frames);
x=vertcat(traj.x);
y=vertcat(traj.y);
z=vertcat(traj.z);
ntraj=vertcat(traj.ntraj);


c=1;
for ktraj=1:numel(traj)
    ntraj(c:c+traj(ktraj).L-1)=traj(ktraj).ntraj;
    c=c+traj(ktraj).L;
end

Nframe=max(frames);

% Sorting of data by their frame number
[frames,ind]=sort(frames);
x=x(ind);
y=y(ind);
z=z(ind);
ntraj=ntraj(ind);


delta_temp=struct([]);
for kf=1:Nframe %for each frame (use parfor to run in parallel)
   % cordinates for frame kf
   index=find(frames==kf);
   xkf=x(index);
   ykf=y(index);
   zkf=z(index);
   ntrajkf=ntraj(index);
   
   % Sort trajectories index
   [ntrajkf,ind] = sort(ntrajkf);
   xkf=xkf(ind);
   ykf=ykf(ind);
   zkf=zkf(ind);
   
   if ~isempty(xkf)
        %separations for frame kf
        L=length(xkf);
        x1=xkf*ones(1,L);
        y1=ykf*ones(1,L);
        z1=zkf*ones(1,L);
        x2=ones(L,1)*xkf';
        y2=ones(L,1)*ykf';
        z2=ones(L,1)*zkf';
        dx=reshape(triu(x2-x1,1),[],1);
        dy=reshape(triu(y2-y1,1),[],1);
        dz=reshape(triu(z2-z1,1),[],1);

        %track index pair for frame kf
        n1=ntrajkf.*ones(1,L);
        n2=ones(L,1).*ntrajkf';
        dn=reshape(triu(n2+1i*n1,1),[],1);

        %remove pairs of same particule
        ind=find(dn==0);
        dx(ind)=[];
        dy(ind)=[];
        dz(ind)=[];
        dn(ind)=[];

        %save data for frame kf
        delta_temp(kf).dx=dx;
        delta_temp(kf).dy=dy;
        delta_temp(kf).dz=dz;
        delta_temp(kf).dt=ones(length(dx),1)*kf;
        delta_temp(kf).dn=dn;
   end
   
    %display frame number every 1000 frames
    if rem(kf,1000)==0
        disp(['Pair dispersion processed for ' num2str(kf) ' frames']);
    end
end

%pair dispersion processing
delta=struct([]);
c=0;

 %sort per track index pair
dx=vertcat(delta_temp.dx);
dy=vertcat(delta_temp.dy);
dz=vertcat(delta_temp.dz);
dt=vertcat(delta_temp.dt);
dn=vertcat(delta_temp.dn);
[dn,ind]=sort(dn,'ComparisonMethod','real');
dx=dx(ind);
dy=dy(ind);
dz=dz(ind);
dt=dt(ind);
        
%rearrange per pair and save data for movie k
pairs=bwconncomp(1-abs(diff(dn))>=1); %identify pairs longer than 2
t0=cellfun(@(x)(x(1)),pairs.PixelIdxList);
Lt=cellfun(@length,pairs.PixelIdxList);
Npair=pairs.NumObjects;
for kp=1:Npair %for each pair
    dxp=dx(t0(kp):t0(kp)+Lt(kp));
    dyp=dy(t0(kp):t0(kp)+Lt(kp));
    dzp=dz(t0(kp):t0(kp)+Lt(kp));
    dp=sqrt(dxp.^2+dyp.^2+dzp.^2);
    delta(kp).dx=dxp;
    delta(kp).dy=dyp;
    delta(kp).dz=dzp;
    delta(kp).d=dp;
    delta(kp).dt=dt(t0(kp):t0(kp)+Lt(kp));
    delta(kp).dn=dn(t0(kp));
    delta(kp).L=length(dp);
end

disp('Pair dispersion processing complete!');
elapsed=toc;
disp(['Elapsed time: ' num2str(elapsed,'%.2f') 's']);
end
