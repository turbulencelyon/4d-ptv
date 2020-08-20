function delta=pairdisp_proc(tracks)

% delta=pairdisp_proc(tracks)
%
% 05/2020 - David Dumont
% 02/2020 - Thomas Basset (adapted from Micka�l Bourgoin)
%
% Pair dispersion processing for positions and velocities
%
% The loop in line 98 can run in parallel to improve time execution.
% _______________________________________________________________
% INPUT
% tracks : track structure with fields {movie,x,y,z,vx,vy,vz,t,L}
%          (field 'movie' required only if more than one movie)
%
% OUTPUTS
% delta(k).dx  : x separation
% delta(k).dy  : y separation
% delta(k).dz  : z separation
% delta(k).d   : total separation
% delta(k).dvx : x velocity difference
% delta(k).dvy : y velocity difference
% delta(k).dvz : z velocity difference
% delta(k).dvp : velocity difference projected along d
% delta(k).dt  : time
% delta(k).dn  : track index pair
% _______________________________________________________________

tic

%coordinates
x=vertcat(tracks.x);
y=vertcat(tracks.y);
z=vertcat(tracks.z);
t=vertcat(tracks.frames);

% %movie and track indices
% n=zeros(length(x),1);
% c=1;
% fields=fieldnames(tracks);
% if ismember('movie',fields)
%     movie=zeros(length(x),1);
%     for k=1:length(tracks)
%         temp=tracks(k).L;
%         movie(c:c+temp-1)=tracks(k).movie;
%         n(c:c+temp-1)=k;
%         c=c+temp;
%     end
% else
%     movie=ones(length(x),1);
%     for k=1:length(tracks)
%         temp=tracks(k).L;
%         n(c:c+temp-1)=k;
%         c=c+temp;
%     end
% end

% %rearrange per movie
% Nmovie=max(movie);
% x=accumarray(movie,x,[Nmovie 1],@(x) {x});
% y=accumarray(movie,y,[Nmovie 1],@(x) {x});
% z=accumarray(movie,z,[Nmovie 1],@(x) {x});
% vx=accumarray(movie,vx,[Nmovie 1],@(x) {x});
% vy=accumarray(movie,vy,[Nmovie 1],@(x) {x});
% vz=accumarray(movie,vz,[Nmovie 1],@(x) {x});
% t=accumarray(movie,t,[Nmovie 1],@(x) {x});
% n=accumarray(movie,n,[Nmovie 1],@(x) {x});

xk = vertcat(tracks.x);
yk = vertcat(tracks.y);
zk = vertcat(tracks.z);
tk = vertcat(tracks.x);

%pair dispersion processing
delta=struct([]);
c=0;
for k=1:Nmovie %for each movie
    
    %coordinates for movie k
    xk=x{k};
    yk=y{k};
    zk=z{k};
    vxk=vx{k};
    vyk=vy{k};
    vzk=vz{k};
    tk=t{k};
    nk=n{k};
    
    if ~isempty(xk)
        %rearrange per frame
        Nframe=max(tk);
        xk=accumarray(tk,xk,[Nframe 1],@(x) {x});
        yk=accumarray(tk,yk,[Nframe 1],@(x) {x});
        zk=accumarray(tk,zk,[Nframe 1],@(x) {x});
        vxk=accumarray(tk,vxk,[Nframe 1],@(x) {x});
        vyk=accumarray(tk,vyk,[Nframe 1],@(x) {x});
        vzk=accumarray(tk,vzk,[Nframe 1],@(x) {x});
        nk=accumarray(tk,nk,[Nframe 1],@(x) {x});
        
        delta_temp=struct([]);
        for kf=1:Nframe %for each frame (use parfor to run in parallel)
            
            %coordinates for frame kf
            xkf=xk{kf};
            ykf=yk{kf};
            zkf=zk{kf};
            vxkf=vxk{kf};
            vykf=vyk{kf};
            vzkf=vzk{kf};
            nkf=nk{kf};
            
            %sort track indices (can be mixed up by previous accumarray)
            [nkf,ind]=sort(nkf);
            xkf=xkf(ind);
            ykf=ykf(ind);
            zkf=zkf(ind);
            vxkf=vxkf(ind);
            vykf=vykf(ind);
            vzkf=vzkf(ind);
            
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
                
%                 %velocity differences for frame kf
%                 vx1=vxkf*ones(1,L);
%                 vy1=vykf*ones(1,L);
%                 vz1=vzkf*ones(1,L);
%                 vx2=ones(L,1)*vxkf';
%                 vy2=ones(L,1)*vykf';
%                 vz2=ones(L,1)*vzkf';
%                 dvx=reshape(triu(vx2-vx1,1),[],1);
%                 dvy=reshape(triu(vy2-vy1,1),[],1);
%                 dvz=reshape(triu(vz2-vz1,1),[],1);
                
                %track index pair for frame kf
                n1=nkf*ones(1,L);
                n2=ones(L,1)*nkf';
                dn=reshape(triu(n2+1i*n1,1),[],1);
                
                %remove pairs of same particule
                ind=find(dn==0);
                dx(ind)=[];
                dy(ind)=[];
                dz(ind)=[];
                dvx(ind)=[];
                dvy(ind)=[];
                dvz(ind)=[];
                dn(ind)=[];
                
                %save data for frame kf
                delta_temp(kf).dx=dx;
                delta_temp(kf).dy=dy;
                delta_temp(kf).dz=dz;
                delta_temp(kf).dvx=dvx;
                delta_temp(kf).dvy=dvy;
                delta_temp(kf).dvz=dvz;
                delta_temp(kf).dt=ones(length(dx),1)*kf;
                delta_temp(kf).dn=dn;
            end
            
            %display frame number every 1000 frames
            if rem(kf,1000)==0
                disp(['Pair dispersion processed for ' num2str(kf) ' frames']);
            end
            
        end
        
        %sort per track index pair
        dx=vertcat(delta_temp.dx);
        dy=vertcat(delta_temp.dy);
        dz=vertcat(delta_temp.dz);
        dvx=vertcat(delta_temp.dvx);
        dvy=vertcat(delta_temp.dvy);
        dvz=vertcat(delta_temp.dvz);
        dt=vertcat(delta_temp.dt);
        dn=vertcat(delta_temp.dn);
        [dn,ind]=sort(dn,'ComparisonMethod','real');
        dx=dx(ind);
        dy=dy(ind);
        dz=dz(ind);
        dvx=dvx(ind);
        dvy=dvy(ind);
        dvz=dvz(ind);
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
            dvxp=dvx(t0(kp):t0(kp)+Lt(kp));
            dvyp=dvy(t0(kp):t0(kp)+Lt(kp));
            dvzp=dvz(t0(kp):t0(kp)+Lt(kp));
            delta(c+kp).dx=dxp;
            delta(c+kp).dy=dyp;
            delta(c+kp).dz=dzp;
            delta(c+kp).d=dp;
            delta(c+kp).dvx=dvxp;
            delta(c+kp).dvy=dvyp;
            delta(c+kp).dvz=dvzp;
            delta(c+kp).dvp=(dvxp.*dxp+dvyp.*dyp+dvzp.*dzp)./dp;
            delta(c+kp).dt=dt(t0(kp):t0(kp)+Lt(kp));
            delta(c+kp).dn=dn(t0(kp));
        end
        c=c+Npair;
        
        disp([num2str(Npair) ' pairs processed for movie ' num2str(k)]); 
    end
    
end

disp('Pair dispersion processing complete!');
elapsed=toc;
disp(['Elapsed time: ' num2str(elapsed,'%.2f') 's']);

end