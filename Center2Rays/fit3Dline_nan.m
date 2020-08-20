function [xyz0,direction]=fit3Dline_nan(XYZ)
%% [xyz0,direction]=fit3Dline_jv(XYZ)
%
% @MBourgoin 01/2019

I = find(isnan(XYZ(:,1)));
XYZ(I,:)=[];

if size(XYZ,1)>2

xyz0=mean(XYZ);
%xyz0=cell2mat(arrayfun(@(x) mean(x.CCrw),Proj,'UniformOutput',false)); 

A=bsxfun(@minus,XYZ,xyz0); %center the data

% xyz0=XYZ(3,:);
% A= XYZ;

% xyz0=XYZ(plan_centre,:);
% A=bsxfun(@minus,XYZ,xyz0); %center the data

%[U,S,V]=svd(A);
[Uac Sac Vac]=arrayfun(@(kkk) svd(A(:,:,kkk)),[1:size(A,3)],'UniformOutput',false);
Ua=cat(3,Uac{:}); 
Sa=cat(3,Sac{:});
Va=cat(3,Vac{:}); clear Uac Sac Vac;

%direction=cross(V(:,end),V(:,end-1));
dd=arrayfun(@(x) cross(Va(:,end,x),Va(:,end-1,x)),[1:size(Va,3)],'UniformOutput',false);
direction=cat(3,dd{:})';  clear dd;
else
    %xyz0 = [NaN NaN NaN];
    %direction = [NaN NaN NaN];
    xyz0=[];
    direction=[];
end

%line = [xyz0'  direction];
