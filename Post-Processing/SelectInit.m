
function [IndexFinal,D]=SelectInit (InitFrame,track,MaxFrame,seuil1,seuil2)
%% Compute the distance between two track at each time for a given initial separation maximum and a given range of frame
%Input:
%InitFrames= first frames in the range
%track=track structure obtain after tracking3D of after stitching
%MaxFrame=last frame of the range 
%Seuil=maximum initial separation 
%Output:
%IndexFinal=First 2 column are the index of the 2 tracks we will observe,
%3rd column is the initial separation and last 2 are the id of the 2 track
%in the track stucture (usefull if you want to visualized the 2 track with visu_1 function )

%D=distance (norm) between the two track
%Dinit=initial seperation for all the track consider in the range 
%------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Supression des trajectoires trop petites 
range=MaxFrame-InitFrame;
%track([track.L]<range)=[];

%% Recherche des indices des trajectoires correspondant à la plage voulu 
index=find([track.L]>range-1);
Ntrack_init=length(index);


index_bis=[];
for i=1:Ntrack_init
    elm=index(i);
    Frames=track(elm).frames;
    Lframes=length(Frames);
    Llimit=find([Frames]==InitFrame);
    Hlimit=find([Frames]==MaxFrame); 
    if ~isempty(Hlimit)
        if ~isempty(Llimit)
            index_bis=[index_bis,elm];
        end
    end
end



%% Ajustement des trajectoires pour qu'elles commencent toutes au même point

Xtt=[];
Ytt=[];
Ztt=[];
Ntrack=length(index_bis);           %Nombre de trajectoire après avoir gardé celle qui correspondent 
fprintf('%d tracks are in this range %d %d frames \n',Ntrack,InitFrame,MaxFrame)
for i=1:Ntrack
    elm=index_bis(i);
    %visu_1(track,elm);
    %hold on
    X=track(elm).x;
    Y=track(elm).y;
    Z=track(elm).z;
    Frames=track(elm).frames;
    minF=min(Frames);
    ajust=InitFrame-minF;
    for k=1:range
        Xaj(k)=X(ajust+k);
        Yaj(k)=Y(ajust+k);
        Zaj(k)=Z(ajust+k);
    end
    Xtt=[Xtt,Xaj];  %On stotck toutes les traj à la suite puisqu'on connait la taille de range 
    Ytt=[Ytt,Yaj];
    Ztt=[Ztt,Zaj];
end 

%% Calcul des distances initiales
Dinit=[];

for i=1:Ntrack-1
    if i==1
        Xinit1=Xtt(1);
        Yinit1=Ytt(1);
        Zinit1=Ztt(1);
    else
        Xinit1=Xtt(i*range+1);
        Yinit1=Ytt(i*range+1);
        Zinit1=Ztt(i*range+1);
    end
    for j=i:Ntrack-1

        if i~=j;
            if j==1
                Xinit2=Xtt(1);
                Yinit2=Ytt(1);
                Zinit2=Ztt(1);
            else
                Xinit2=Xtt(j*range+1);
                Yinit2=Ytt(j*range+1);
                Zinit2=Ztt(j*range+1);   
            end
            Dinit(i,j)=sqrt((Xinit1-Xinit2)^2+(Yinit1-Yinit2)^2+(Zinit1-Zinit2)^2);
        end
    end
end


%% Recherche des distance initial inférieur à un certain seuil
L=[];
if ~isempty(Dinit)
    [L,C]=find(Dinit<seuil2 & Dinit>seuil1);
    Coord=[C,L];
    IndexFinal=[];
    for i=1:length(L)
        A=Coord(i,1);
        B=Coord(i,2);
        C=[A,B,Dinit(B,A),index_bis(A),index_bis(B)];
        IndexFinal=vertcat(IndexFinal,C);
    end

    fprintf('%d track have an initial separation between %dmm and %dmm for the frames %d to %d \n',length(L),seuil1,seuil2,InitFrame,MaxFrame)

    %% Calcul des incréments de vitesse au cour du temps 
    NtProche=length(L);
    hold on 
    if NtProche ~=0
            for i=1:NtProche
                if IndexFinal(i,1)==1
                    A1=1;
                    Xa=Xtt(A1:A1+range);
                    Ya=Ytt(A1:A1+range);
                    Za=Ztt(A1:A1+range);
                else
                    A1=IndexFinal(i,1)*range;
                    Xa=Xtt(A1+1:A1+range);
                    Ya=Ytt(A1+1:A1+range);
                    Za=Ztt(A1+1:A1+range);
                end
                if IndexFinal(i,2)==1
                    B1=1;
                    Xb=Xtt(B1:B1+range);
                    Yb=Ytt(B1:B1+range);
                    Zb=Ztt(B1:B1+range);
                else
                    B1=IndexFinal(i,2)*range;
                    Xb=Xtt(B1+1:B1+range);
                    Yb=Ytt(B1+1:B1+range);
                    Zb=Ztt(B1+1:B1+range);
                end
            %A1=IndexFinal(i,1)*range
            %B1=IndexFinal(i,2)*range;

                for j=1:range
                    D(i,j)=((Xa(j)-Xb(j))^2+(Ya(j)-Yb(j))^2+(Za(j)-Zb(j))^2);
                end
            end


            for i=1:NtProche
                %plot(D(i,:))
            end
    end
else
    D=[];
    IndefFinal=[];
end
if length(L)==0
    D=[];
    IndexFinal=[];
end
end


 
    