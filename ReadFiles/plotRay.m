function plotRay(P,V,s,varargin)

Nrays = size(P,1);
hold on
for k = 1:Nrays
    X = (P(k,:) + V(k,:).*s');
    X = sort(X,3);
    
    plot3(X(:,1),X(:,2),X(:,3),varargin{:});
end
%hold off

