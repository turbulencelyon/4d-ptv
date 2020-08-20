function S=SelectLocation(track,x1,x2,y1,y2,z1,z2)


L=length(track);
S=[];
index=[];
for i=1:L
    X=track(i).x;
    Y=track(i).y;
    Z=track(i).z;
    A=any(X>x1 & X<x2);
    B=any(Y>y1 & Y<y2);
    C=any(Z>z1 & Z<z2);
    
    if A==1 && B==1 && C==1
        disp(i);
        index=[index,i];
    end
end
S=track(index);
end