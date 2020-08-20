function mcell=meancell(cellarray,varargin)

% mcell=meancell(cellarray,varargin)
%
% 02/2020 - Thomas Basset (adapted from Mickaël Bourgoin)
%
% Mean of array cells of different lengths
% ______________________________________________
% INPUTS
% cellarray : array cells of different lengths
% weights   : weights to compute mean (optional)
%
% OUTPUTS
% mcell(k).mean  : mean
% mcell(k).std   : standard deviation
% mcell(k).Npts  : number of points (if weights)
% mcell(k).Ncell : number of cells
% ______________________________________________

%weights
if nargin>1
    weights=varargin{1};
end

%means
L=cellfun(@length,cellarray);
for k=1:max(L)
    J=find(L>=k);
    if nargin>1
        Npts=nansum(cellfun(@(y)(y(k)),weights(J)));
        mcell.mean(k)=nansum(cellfun(@(x,y)(x(k)*y(k)),cellarray(J),weights(J)))/Npts;
        mcell.std(k)=sqrt(nansum(cellfun(@(x,y)((x(k)-mcell.mean(k))^2*y(k)),cellarray(J),weights(J)))/Npts);
        mcell.Npts(k)=Npts;
    else
        Npts=length(J);
        mcell.mean(k)=nansum(cellfun(@(x)(x(k)),cellarray(J)))/Npts;
        mcell.std(k)=sqrt(nansum(cellfun(@(x)((x(k)-mcell.mean(k))^2),cellarray(J)))/Npts);
    end
    mcell.Ncell(k)=length(J);
end

end