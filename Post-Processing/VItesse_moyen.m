
index=find([tracks.L]);
dt=1/150;
V_tt=0;
for i=1 : length(index)
    elm=index(i);
    X=tracks(elm).x;   
    Y=tracks(elm).y;
    Z=tracks(elm).z;
    Sum_V=0;
    longueur_V=length(X)-1;
    for j=1 : longueur_V
        Vx(j)=(X(j)-X(j+1))/dt;
        Vy(j)=(X(j)-X(j+1))/dt;
        Vz(j)=(X(j)-X(j+1))/dt;
        V(j)=Vx(j)+Vy(j)+Vz(j)his;
        Sum_V=Sum_V+V(j);
    end
    V_mean=Sum_V/longueur_V;
    V_tt(i)=V_mean;
end
%Vseuil=0;
%indexbis=find([V_tt]>5);
%for i=1 : length(indexbis)
%Vseuil(i)=V_tt(indexbis(i));    
    
end


mean(V_tt)
std(V_tt)
