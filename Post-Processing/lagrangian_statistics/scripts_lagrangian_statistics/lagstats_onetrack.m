function [r,s,Npts]=lagstats_onetrack(track,n,varargin)

% [r,s,Npts]=lagstats_onetrack(track,n,varargin)
%
% 02/2020 - Thomas Basset (adapted from Mickaël Bourgoin)
%
% Lagrangian statistics for one track (autocorrelation and structure functions)
% _____________________________________________________
% INPUTS
% track : track array
% n     : structure function order
% t     : time array (optional, required if time holes)
%
% OUTPUTS
% r    : autocorrelation function
% s    : n order structure function
% Npts : number of points
% _____________________________________________________

%regular time increments
if nargin>2
    t=varargin{1};
    t=t-t(1)+1;
    if t(end)~=length(t) %if time holes
        temp=NaN(t(end),1);
        temp(t)=track;
        track=temp;
    end
end
t=(1:length(track))-1;

%compute functions
r_cell=cell(length(t),1);
s_cell=cell(length(t),1);
for k=1:length(t) %for each time increment
    ind=(1:(length(track)-t(k)));
    r_cell{k}=track(ind+t(k)).*track(ind); %autocorrelation function
    s_cell{k}=(track(ind+t(k))-track(ind)).^n; %n order structure function
end
r=cellfun(@nanmean,r_cell);
s=cellfun(@nanmean,s_cell);
Npts=cellfun(@(x)(length(x)-sum(isnan(x))),r_cell);

end