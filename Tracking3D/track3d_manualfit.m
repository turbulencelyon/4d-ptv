function [tracks,traj]=track3d_manualfit(folderout,FileName,data,maxdist,lmin,flag_pred,npriormax,flag_conf)

% 03/2019 - Thomas Basset 
% 01/2020 - "Optimisation" David Dumont
%
% Track particles on 3D images
% ____________________________________________________________________________
% INPUTS
% folderout : path to folder where track results will be saved
% FileName  : Rays file name 
% data      : array composed of all matches [FrameNumber, x, y, z, Error]
% maxdist   : maximum travelled distance between two successive frames
% lmin      : minimum length of a trajectory (number of frames)
% flag_pred : 1 for predictive tracking, 0 otherwise
% npriormax : maximum number of prior frames used for predictive tracking
% flag_conf : 1 for conflict solving, 0 otherwise
%
% OUTPUTS
% traj(kt).ntraj  : trajectory index
% traj(kt).L      : trajectory length
% traj(kt).frames : trajectory frames
% traj(kt).x      : x-position
% traj(kt).y      : y-position
% traj(kt).z      : z-position
% traj(kt).nmatch : element indices in tracks
% tracks          : trajectory raw data
% ____________________________________________________________________________

fprintf('Beginning tracking...\n');

%raw data
tracks=zeros(size(data,1),7);
tracks(:,1:4)=data(:,1:4); %frame number and position
tracks(:,5)=0; %trajectory index
tracks(:,6)=1; %availability of particles (0 unavailable, 1 available, 2 conflict)
tracks(:,7)=(1:length(data(:,1)))'; %element index

%number of particles in first frame = number of active trajectories
ind_act=find(tracks(:,1)==min(tracks(:,1)));
tracks(ind_act,5)=ind_act;
tracks(ind_act,6)=0;
ntraj=max(ind_act);

% Array containing in the line n the index of particules which belong to
% trajectory n.
partintraj = struct();
for i= 1:max(ind_act)
    partintraj(i).partindex(1) = i;
end

dispmax=1e5;
tic
fprintf("Let's start loop over frame...\n")
for kf=min(tracks(:,1))+1:max(tracks(:,1)) %for each frame
    
    ind_new=find(tracks(:,1)==kf);
    min_new=min(ind_new);
    max_new=max(ind_new);
    c=1;
    dispoind=zeros(dispmax,1);
    dispoact=zeros(dispmax,1);
    disposqrd=zeros(dispmax,1);
    
    %positions of particles in frame kf
    newx=tracks(ind_new,2);
    newy=tracks(ind_new,3);
    newz=tracks(ind_new,4);

    for kp=1:length(ind_act) %for each particle in frame kf-1
        
        %particle kp in frame kf-1
        actx=tracks(ind_act(kp),2); %x-position
        acty=tracks(ind_act(kp),3); %y-position
        actz=tracks(ind_act(kp),4); %z-position 
        actnum=tracks(ind_act(kp),5); %trajectory index
        
        %predictive tracking
        if flag_pred==1
            %indices of up to last npriormax points in active trajectory
            if numel(partintraj(actnum).partindex)<npriormax
                ind_pred = partintraj(actnum).partindex;
            else
                ind_pred = partintraj(actnum).partindex(end-npriormax+1:end);
            end
            %length of used active trajectory part
            nprior=length(ind_pred);
            if nprior>1
                xmin = mean(tracks(ind_pred(1:floor(nprior/2)),2));
                ymin = mean(tracks(ind_pred(1:floor(nprior/2)),3));
                zmin = mean(tracks(ind_pred(1:floor(nprior/2)),4));
                
                xmax = mean(tracks(ind_pred(end-floor(nprior/2)+1:end),2));
                ymax = mean(tracks(ind_pred(end-floor(nprior/2)+1:end),3));
                zmax = mean(tracks(ind_pred(end-floor(nprior/2)+1:end),4));
                
                actx = xmax+(xmax-xmin)/(floor((nprior-1)/2)+1);
                acty = ymax+(ymax-ymin)/(floor((nprior-1)/2)+1);
                actz = zmax+(zmax-zmin)/(floor((nprior-1)/2)+1);
            end
        end
        
        %calculation of squared distance: nearest neighbour
        sqrdist=(actx-newx).^2+(acty-newy).^2+(actz-newz).^2;
        
        %minimize the distance
        [sqrdmin,ind_min]=min(sqrdist);
        if sqrdmin<maxdist^2
            dispo=tracks(ind_new(ind_min),6);
            
            %available particle
            if dispo==1 
                %particle becomes unavailable
                tracks(ind_new(ind_min),6)=0;
                %new trajectory index
                tracks(ind_new(ind_min),5)=actnum;
                partintraj(actnum).partindex(end+1) = ind_new(ind_min);
                %useful arrays for next steps
                dispoind(c)=ind_new(ind_min); %element index
                dispoact(c)=actnum; %trajectory index
                disposqrd(c)=sqrdmin; %minimum squared distance
                c=c+1;
                
            %unavailable particle = conflict
            elseif dispo==0
                i=find(dispoind==ind_new(ind_min));
                actnump=dispoact(i);
                
                %non solving
                if (flag_conf==0 || sqrdmin==disposqrd(i))
                    %particle could go to two different trajectories
                    tracks(ind_new(ind_min),6)=2;
                    %trajectory index set to zero 
                    tracks(ind_new(ind_min),5)=0;
                    partintraj(actnump).partindex(end) = [];
                end
                
                %solving
                if (flag_conf==1 && sqrdmin<disposqrd(i)) %new track is a better match
                    %new trajectory index
                    tracks(ind_new(ind_min),5)=actnum;
                    partintraj(actnum).partindex(end+1) = ind_new(ind_min);
                    partintraj(actnump).partindex(end) = [];
                    %array updating
                    dispoact(i)=actnum;
                    disposqrd(i)=sqrdmin;
                end

            end
            
        end
        
    end
    
    %tracked particules: found once or not found
    tracks_temp=tracks(min_new:max_new,:);
    %found once
    ind=tracks_temp(:,6)==0;
    ind_act0=tracks_temp(ind,7);
    %not found
    ind=tracks_temp(:,6)==1;
    ind_act1=tracks_temp(ind,7);
    if ~isempty(ind_act1) %new trajectory index for particles not found
        for k=1:length(ind_act1)
            ntraj=ntraj+1;
            tracks(ind_act1(k),5)=ntraj;
            partintraj(ntraj).partindex(1) = ind_act1(k);
        end
    end
    ind_act=[ind_act0;ind_act1];
    
    %display frame number every 100 frames
    if rem(kf,100)==0 
        disp(['Trajectories found for ' num2str(kf) ' frames']);
    end
    
end

fprintf("Tracking finished. Structure creation in progress...")

%% structure creation
tracks_array = sortrows(tracks,5);
% remove trajecti
fstart = zeros(length(tracks_array),1);
fstart(1) = 1; fstart(2:end) = diff(tracks_array(:,5));
fstart = find(fstart==1); % index where new trajectory starts

flong = diff(fstart); % length of each trajectory
last_traj = find(tracks_array(:,5)==length(fstart));
flong(end+1) = length(last_traj);

ftrackselec = find(flong>=lmin); % index of trajectories longer than lmin
traj = struct('ntraj',cell(1,numel(i)),'L',[],'frames',[],'x',[],'y',[],'z',[],'nmatch',[]);
for kk = 1 :length(ftrackselec)
    ll  = flong(ftrackselec(kk));
    deb = fstart(ftrackselec(kk));
    traj(kk).L = ll;
    traj(kk).ntraj = kk;
    traj(kk).frames = tracks_array(deb:deb+ll-1,1);
    traj(kk).nmatch = tracks_array(deb:deb+ll-1,7);  % ki previously
    traj(kk).x = tracks_array(deb:deb+ll-1,2);
    traj(kk).y = tracks_array(deb:deb+ll-1,3);
    traj(kk).z = tracks_array(deb:deb+ll-1,4);
end

disp([num2str(kk) ' trajectories longer than ' num2str(lmin) ' frames (from ' num2str(length(vertcat(traj.x))) ' matches)']);
fprintf('Saving to .h5 file in %s\n', folderout);

%save in .h5
FileSaveName = fullfile(folderout,['tracks_' char(FileName) '.h5']);
fields=fieldnames(traj);
for k=1:length(fields) %for each field
    fieldk=fields{k};
    h5create(FileSaveName,['/' fieldk],[length(vertcat(traj.(fieldk))) 1],'ChunkSize',[1000 1],'Deflate',9)
    h5write(FileSaveName,['/' fieldk],vertcat(traj.(fieldk)))
end
% saving of tracking parameters
for param=["maxdist","lmin","flag_pred","npriormax","flag_conf"]
    h5create(FileSaveName,['/trackingparameters/' char(param)],1)
    h5write(FileSaveName,['/trackingparameters/' char(param)],eval(param))
end


disp('Tracking complete!');
elapsed=toc;
disp(['Elapsed time: ',num2str(elapsed,'%.2f'),'s']);

end
