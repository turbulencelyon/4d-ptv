function [Struct,XAx]=StructureFct(SeuilInit,SeuilFinal,DeltaSeuil,Pas,NbFrames,Epsilon,Ech,Dinit,S,index,Comp)
Time=linspace(0,NbFrames/Ech,NbFrames);
Delta=[];
Struct=[];
L=1;
Ltot=[];
figure

for i=SeuilInit:Pas:SeuilFinal
    A=find(Dinit<i+DeltaSeuil & Dinit>i);
    L=length(A)
    while L<1500
        DeltaSeuil=DeltaSeuil+0.002
        A=find(Dinit<i+DeltaSeuil & Dinit>i);
        L=length(A);
    end
    L=1;
    Tstar=((((2*i+DeltaSeuil)/2)*10^-3)^2/(Epsilon))^(1/3);
    TimeNormalized=Time./Tstar; 
    Lim=length(find(TimeNormalized<0.1));
    
    
%% Call the function that will compute the distance between two tracks
    Delta=SelectSizeComp(S,Dinit,index,i,i+DeltaSeuil,Comp);
    Dm=mean(Delta)*10^-6;   % conversion mm2 => m^2
    TimeNorm=Time.^2;
    Dmnorm=Dm(1:1000)./TimeNorm;
    loglog(Time(20:Lim),Dmnorm(20:Lim),'-O','DisplayName',sprintf('Initial separation %.2e-%.2e mm',i,(i+DeltaSeuil)));
    fo=fitoptions('Method','NonLinearLeastSquare','Lower',[0],'Upper',[max(Dm)],'StartPoint',[0]);
    ft=fittype('a*x^n+c','problem','n');
    fitting=fit(Time(20:Lim)',Dm(20:Lim)',ft,'problem',2);
    %[X,s]=polyfit(Time(5:Lim),Dm(5:Lim),2);
    X=coeffvalues(fitting)
    loglog(Time(1:Lim),X(1)*Time(1:Lim).^2,'x')
    hold on; 
    XAx(i)=i+DeltaSeuil/2;

    Struct=[Struct,X(1)];
end 

XAx(XAx==0)=[];

%% Graph parameters
%xlabel('t/t*');
%ylabel('<D_0^2>(mm²)/t^{3}\epsilon');
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
%legend('location','best');
%grid on;

end