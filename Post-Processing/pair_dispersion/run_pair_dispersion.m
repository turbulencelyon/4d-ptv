%% Pair dispersion processing
delta=pairdisp_proc(tracks);



%% Pair dispersion statistics for separations
bin_type='log';
dmin=1;
dmax=50;
Nbin=20;
delta_d=pairdisp_stats_d(delta,bin_type,dmin,dmax,Nbin);

%% Time evolution
fps=6250;
d0=vertcat(delta_d.d0);
figure;hold on;grid on
for k=1:length(delta_d)
    y=delta_d(k).d2; %mean square separation
    %y=delta_d(k).Npts; %number of pairs
    x=(0:length(y)-1)'/fps; %time
    %plot(y,'x')
    %plot(x,y,'x')
    plot(y./x.^2,'x') %compensated
    %plot(x./d0(k)^(2/3),y./d0(k)^2,'x') %normalized
end
set(gca,'XScale','log','YScale','log')
legend(num2str(d0),'Location','northwest')

%% Structure function
S2=zeros(length(delta_d),1);
for k=1:length(S2)
    y=delta_d(k).d2;
    x=(0:length(y)-1)'/fps;
    y_comp=y./x.^2;
    S2(k)=median(y_comp(10:50))*1e-6; %from mm to m
end
d0_m=d0*1e-3; %from mm to m
figure;loglog(d0_m,S2,'x');grid on

%% Epsilon
epsilon=((3/22)*S2.^1.5).*d0_m.^-1;
figure;loglog(d0_m,epsilon,'x');grid on



%% Pair dispersion statistics for velocities
bin_type='lin';
dmin=10;
dmax=30;
Nbin=10;
delta_v=pairdisp_stats_v(delta,bin_type,dmin,dmax,Nbin);

%% Time evolution
fps=6250;
d0=vertcat(delta_v.d0);
figure;hold on;grid on
for k=1:length(delta_v)
    y=(delta_v(k).Rx+delta_v(k).Ry+delta_v(k).Rz)/3; %autocorrelation function
    %y=delta_v(k).Npts; %number of points
    x=(0:length(y)-1)'/fps; %time
    y=y(1:350);
    x=x(1:350);
    %plot(y,'x')
    %plot(x,y,'x')
    plot(y/y(1),'x')
    %plot(x/(d0(k)^(2/3)),y/y(1),'x') %normalized
end
set(gca,'YScale','log')
legend(num2str(d0),'Location','northeast')

%% Time scale: exponential fit
coeff=zeros(length(delta_v),2);
figure;hold on;grid on
for k=1:length(delta_v)
    y=(delta_v(k).Rx+delta_v(k).Ry+delta_v(k).Rz)/3;
    x=(0:length(y)-1)'./fps;
    y=y(50:150);
    x=x(50:150);
    xyfit=polyfit(x,log(y),1);
    a=xyfit(1);
    b=xyfit(2);
    coeff(k,:)=[a b];
    plot(x,y,'x')
    plot(x,exp(a*x+b),'-');
end
T=-coeff(1:end,1).^-1;

%% Time scale: integral
T=zeros(length(delta_v),1);
figure;hold on;grid on
for k=1:length(delta_v)
    y=(delta_v(k).Rx+delta_v(k).Ry+delta_v(k).Rz)/3;
    x=(0:length(y)-1)'./fps;
    y=y/y(1);
    y=y(1:350);
    x=x(1:350);
    plot(x,y,'x')
    xyfit=polyfit(x,y,5); 
    fun=@(t)xyfit(1)*t.^5+xyfit(2)*t.^4+xyfit(3)*t.^3+xyfit(4)*t.^2+xyfit(5)*t+xyfit(6);
    T(k)=integral(fun,0,x(end));
    plot(x,polyval(xyfit,x),'-')
end

%% Fit T
x=d0(1:end);
y=T(1:end);
xyfit=polyfit(log(x),log(y),1);
a=xyfit(1);
b=xyfit(2);
figure;plot(x,y,'-x');grid on
hold on;plot(x,exp(b)*x.^a,'-')

%% Alpha
epsilon=0.9;
d0_m=d0*1e-3; %from mm to m
T0=(22/6)*epsilon^(-1/3)*d0_m.^(2/3);
alpha=T./T0;

%% Fit alpha
x=d0(1:end);
y=alpha(1:end);
xyfit=polyfit(log(x),log(y),1);
a=xyfit(1);
b=xyfit(2);
figure;plot(x,y,'-x');grid on
hold on;plot(x,exp(b)*x.^a,'-')
