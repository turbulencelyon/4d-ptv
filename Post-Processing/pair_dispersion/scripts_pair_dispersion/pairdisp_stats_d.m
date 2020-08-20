function delta_d=pairdisp_stats_d(delta,bin_type,dmin,dmax,Nbin)

% delta_d=pairdisp_stats_d(delta,bin_type,dmin,dmax,Nbin)
%
% 02/2020 - Thomas Basset (adapted from Mickaël Bourgoin)
%
% Pair dispersion statistics for separations
% _________________________________________________________
% INPUTS
% delta    : structure returned by "pairdisp_proc"
% bin_type : 'lin' ('log') for linear (logarithmic) binning
% dmin     : minimum initial separation
% dmax     : maximum initial separation
% Nbin     : number of bins for initial separation
%
% OUTPUTS
% delta_d(k).dx2  : mean square x separation
% delta_d(k).dy2  : mean square y separation
% delta_d(k).dz2  : mean square z separation
% delta_d(k).d2   : mean square total separation
% delta_d(k).Npts : number of pairs
% delta_d(k).d0   : initial total separation
% _________________________________________________________

tic

%initial separation histogram
%dx0=arrayfun(@(x)(x.dx(1)),delta);
%dy0=arrayfun(@(x)(x.dy(1)),delta);
%dz0=arrayfun(@(x)(x.dz(1)),delta);
d0=arrayfun(@(x)(x.d(1)),delta);
if strcmp(bin_type,'lin')
    edged=linspace(dmin,dmax,Nbin+1); 
end
if strcmp(bin_type,'log')
    edged=logspace(log10(dmin),log10(dmax),Nbin+1);
end
%[Ndx0,edgesdx0,bindx0]=histcounts(dx0,edged);
%[Ndy0,edgesdy0,bindy0]=histcounts(dy0,edged);
%[Ndz0,edgesdz0,bindz0]=histcounts(dz0,edged);
%[Nd0,edgesd0,bind0]=histcounts(d0,edged);
[~,~,bind0]=histcounts(d0,edged);

%pair dispersion statistics
delta_d=struct([]);
edged0=(edged(1:end-1)+edged(2:end))/2;
for k=1:Nbin %for each initial separation
    
    %substract initial separation
    delta_temp=delta(bind0==k);
    dx2=arrayfun(@(x)((x.dx-x.dx(1)).^2),delta_temp,'UniformOutput',false);
    dy2=arrayfun(@(x)((x.dy-x.dy(1)).^2),delta_temp,'UniformOutput',false);
    dz2=arrayfun(@(x)((x.dz-x.dz(1)).^2),delta_temp,'UniformOutput',false);
    d2=arrayfun(@(x)((x.dx-x.dx(1)).^2+(x.dy-x.dy(1)).^2+(x.dz-x.dz(1)).^2),delta_temp,'UniformOutput',false);
    dt=arrayfun(@(x)(x.dt),delta_temp,'UniformOutput',false);
    
    %same lengths to compute means
    Ld=arrayfun(@(x)(x.dt(end)-x.dt(1)+1),delta_temp);
    Ldmax=max(Ld);
    Npair=length(Ld);
    ndx2=NaN(Npair,Ldmax);
    ndy2=NaN(Npair,Ldmax);
    ndz2=NaN(Npair,Ldmax);
    nd2=NaN(Npair,Ldmax);
    Npts=NaN(Npair,Ldmax);
    for kp=1:Npair
        t=dt{kp}; %consider time holes
        t=t-t(1)+1;
        ndx2(kp,t)=dx2{kp};
        ndy2(kp,t)=dy2{kp};
        ndz2(kp,t)=dz2{kp};
        nd2(kp,t)=d2{kp};
        Npts(kp,t)=ones(1,length(t));
    end
    
    %means
    dx2=nanmean(ndx2);
    dy2=nanmean(ndy2);
    dz2=nanmean(ndz2);
    d2=nanmean(nd2);
    Npts=nansum(Npts);
    
    %save data
    delta_d(k).dx2=dx2';
    delta_d(k).dy2=dy2';
    delta_d(k).dz2=dz2';
    delta_d(k).d2=d2';
    delta_d(k).Npts=Npts';
    delta_d(k).d0=edged0(k);
    
    disp(['Statistics done for bin ' num2str(k)]);
    
end

disp('Pair dispersion statistics complete!');
elapsed=toc;
disp(['Elapsed time: ' num2str(elapsed,'%.2f') 's']);

end