%% Lagrangian statistics
n=2;
[R,S,Npts,Ntrack]=lagstats_tracks(tracks,'vx',n,'t');



%% Autocorrelation function
fps=6250;
t=(0:length(R)-1)'/fps;
figure;plot(t,R,'x');grid on
set(gca,'YScale','log')

%% Time scale: exponential fit
x=t(100:200);
y=R(100:200);
xyfit=polyfit(x,log(y),1);
a=xyfit(1);
b=xyfit(2);
T=-1/a;
figure;plot(y,'x');grid on
hold on;plot(exp(a*x+b),'-')

%% Time scale: integral
x=t(1:400);
y=R(1:400);
y=y/y(1);
xyfit=polyfit(x,y,5);
fun=@(t)(xyfit(1)*t.^5+xyfit(2)*t.^4+xyfit(3)*t.^3+xyfit(4)*t.^2+xyfit(5)*t+xyfit(6));
T=integral(fun,0,x(end));
figure;plot(x,y,'x');grid on
hold on;plot(x,polyval(xyfit,x),'-')



%% Second order structure function
fps=6250;
t=(0:length(S)-1)'/fps;
S=S*(1e-3)^n; %from mm to m
figure;plot(t,S,'x');grid on
set(gca,'XScale','log','YScale','log')

%% Compensated inertial range
figure;plot(t,S./t,'x');grid on
set(gca,'XScale','log','YScale','log')
