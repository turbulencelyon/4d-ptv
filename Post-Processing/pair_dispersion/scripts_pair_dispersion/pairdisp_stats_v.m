function delta_v=pairdisp_stats_v(delta,bin_type,dmin,dmax,Nbin)

% delta_v=pairdisp_stats_v(delta,bin_type,dmin,dmax,Nbin)
%
% 02/2020 - Thomas Basset (adapted from Mickaël Bourgoin)
%
% Pair dispersion statistics for velocities
% _________________________________________________________
% INPUTS
% delta    : structure returned by "pairdisp_proc"
% bin_type : 'lin' ('log') for linear (logarithmic) binning
% dmin     : minimum initial separation
% dmax     : maximum initial separation
% Nbin     : number of bins for initial separation
%
% OUTPUTS
% delta_v(k).Rx    : dvx autocorrelation function
% delta_v(k).Ry    : dvy autocorrelation function
% delta_v(k).Rz    : dvz autocorrelation function
% delta_v(k).Rp    : dvp autocorrelation function
% delta_v(k).Sx    : dvx second order structure function
% delta_v(k).Sy    : dvy second order structure function
% delta_v(k).Sz    : dvz second order structure function
% delta_v(k).Sp    : dvp second order structure function
% delta_v(k).Npts  : number of points
% delta_v(k).Npair : number of pairs
% delta_v(k).d0    : initial total separation
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
delta_v=struct([]);
edged0=(edged(1:end-1)+edged(2:end))/2;
for k=1:Nbin %for each initial separation
    
    %compute functions
    delta_temp=delta(bind0==k);
    [Rx,Sx,Npts,Npair]=lagstats_tracks(delta_temp,'dvx',2,'dt');
    [Ry,Sy,~,~]=lagstats_tracks(delta_temp,'dvy',2,'dt');
    [Rz,Sz,~,~]=lagstats_tracks(delta_temp,'dvz',2,'dt');
    [Rp,Sp,~,~]=lagstats_tracks(delta_temp,'dvp',2,'dt');
    
    %save data
    delta_v(k).Rx=Rx;
    delta_v(k).Ry=Ry;
    delta_v(k).Rz=Rz;
    delta_v(k).Rp=Rp;
    delta_v(k).Sx=Sx;
    delta_v(k).Sy=Sy;
    delta_v(k).Sz=Sz;
    delta_v(k).Sp=Sp;
    delta_v(k).Npts=Npts;
    delta_v(k).Npair=Npair;
    delta_v(k).d0=edged0(k);
    
    disp(['Statistics done for bin ' num2str(k)]);
    
end

disp('Pair dispersion statistics complete!');
elapsed=toc;
disp(['Elapsed time: ' num2str(elapsed,'%.2f') 's']);

end