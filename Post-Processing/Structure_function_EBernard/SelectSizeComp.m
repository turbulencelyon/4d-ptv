
function [Delta]=SelectSizeComp(S,Dinit,index,seuil1,seuil2,Comp)

%% Compute the distance between two track for all tracks that have the rigth initial separation  

%% The input S,Dinit,index come from the function Dinitial.m
% seuil1        : the minimum separation initial to be count
% seuil2        : the maximum separation initial to be count 
%
%% Output
% delta         : a matrix containing the pair dispersion for every couple
% of tracks 
%%

%% Find the track that correspond to the initial separation we want 
L=[];
if ~isempty(Dinit)
    [L,C]=find(Dinit<seuil2 & Dinit>seuil1);
    Coord=[C,L];
    IndexFinal=[];
    for i=1:length(L)
        A=Coord(i,1);
        B=Coord(i,2);
        C=[A,B,Dinit(B,A),index(A),index(B)];
        IndexFinal=vertcat(IndexFinal,C);
    end

    fprintf('%d track have an initial separation between %dmm and %dmm \n',length(L),seuil1,seuil2)

 


%% Compute the separation for each couple of particule in function of time 
    NtProche=length(L);
    hold on 
    if NtProche ~=0
        for i=1:NtProche
            if rem(i,1000)==0
                disp(sprintf('%d %s done',i/NtProche*100,'%'));
            end
            if Comp=='x'
                Xa=S(IndexFinal(i,1)).X;
                Xb=S(IndexFinal(i,2)).X;
                LimA=length(Xa);
                LimB=length(Xb);
                Min=min(LimA,LimB);
                Xa=Xa(1:Min);
                Xb=Xb(1:Min);
                D(i,1:Min)=(Xa-Xb);       % Delta_k(t)
                Delta(i,1:Min)=(D(i,1:Min)-D(i,1)).^2;      % (Delta_k(t)-delta_k(0))^2
            end
            if Comp=='y'
                Ya=S(IndexFinal(i,1)).Y;
                Yb=S(IndexFinal(i,2)).Y;
                LimA=length(Ya);
                LimB=length(Yb);
                Min=min(LimA,LimB);
                Ya=Ya(1:Min);
                Yb=Yb(1:Min);
                D(i,1:Min)=(Ya-Yb);
                Delta(i,1:Min)=(D(i,1:Min)-D(i,1)).^2;
            end
            if Comp=='z'
                Za=S(IndexFinal(i,1)).Z;
                Zb=S(IndexFinal(i,2)).Z;
                LimA=length(Za);
                LimB=length(Zb);
                Min=min(LimA,LimB);
                Za=Za(1:Min);
                Zb=Zb(1:Min);
                D(i,1:Min)=(Za-Zb);
                Delta(i,1:Min)=(D(i,1:Min)-D(i,1)).^2;
            end
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


 
    