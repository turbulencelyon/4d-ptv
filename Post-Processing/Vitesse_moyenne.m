function X=Vitesse_moyenne(tracks,range) 
index=find([tracks.L]);
dt=1/150;
V_tt=0;
for i=1 : length(index)
    elm=index(i);
    X=tracks(elm).x;   
    Y=tracks(elm).y;
    Z=tracks(elm).z;
    Sum_Vx=0;
    Sum_Vy=0;
    Sum_Vz=0;
    longueur_V=length(X)-1;
    for j=1 : longueur_V
        Vx(j)=(X(j)-X(j+1))/dt;
        Vy(j)=(Y(j)-Y(j+1))/dt;
        Vz(j)=(Z(j)-Z(j+1))/dt;
        Sum_Vx=Sum_Vx+Vx(j);
        Sum_Vy=Sum_Vy+Vy(j);
        Sum_Vz=Sum_Vz+Vz(j);
        
    end
    V_meanx=Sum_Vx/longueur_V;
    V_ttx(i)=V_meanx;
    V_meany=Sum_Vy/longueur_V;
    V_tty(i)=V_meany;
    V_meanz=Sum_Vz/longueur_V;
    V_ttz(i)=V_meanz;
end 

%histogram(V_ttx,500)
 hold on 
 
X=histogram(V_tty,500)
%X=histogram(V_ttz,range)
legend('X')
xlabel('Vitesse')
ylabel('Nb traj')
end 
