function [R,S,Npts,Ntrack]=lagstats_tracks(tracks,fieldv,n,varargin)

% [R,S,Npts,Ntrack]=lagstats_tracks(tracks,fieldv,n,varargin)
%
% 02/2020 - Thomas Basset (adapted from Mickaël Bourgoin)
%
% Lagrangian statistics for a track set (autocorrelation and structure functions)
% ___________________________________________________________
% INPUTS
% tracks : track structure
% fieldv : velocity field name
% n      : structure function order
% fieldt : time field name (optional, required if time holes)
%
% OUTPUTS
% R      : autocorrelation function
% S      : n order structure function
% Npts   : number of points
% Ntrack : number of tracks
% ___________________________________________________________

%autocorrelation and n order structure function for each track
if nargin>3
    fieldt=varargin{1};
    [r,s,Npts]=arrayfun(@(x)(lagstats_onetrack(x.(fieldv),n,x.(fieldt))),tracks,'UniformOutput',false);
else
    [r,s,Npts]=arrayfun(@(x)(lagstats_onetrack(x.(fieldv),n)),tracks,'UniformOutput',false);
end

%mean over all tracks
R=meancell(r,Npts);
S=meancell(s,Npts);
Npts=(R.Npts)';
Ntrack=(R.Ncell)';
R=(R.mean)';
S=(S.mean)';

end