%% https://www.mathworks.com/help/signal/ug/hilbert-transform.html
clear; close all;
t = 0:1/1024:1;
x = sin(2*pi*60*t);

%********
figure;
plot(t(1:50),x(1:50),'r');hold on;

y = hilbert(x);

plot(t(1:50),real(y(1:50)))
hold on
plot(t(1:50),imag(y(1:50)))
hold off
axis([0 0.05 -1.1 2])
legend('Real Part','Imaginary Part')


%******
y2=hilbert(x);
plot(t,y);
hold on;
plot(t,real(y),'g');
plot(t,y2,'m');
%% https://www.mathworks.com/help/signal/ug/envelope-extraction-using-the-analytic-signal.html
clear;

t = 0:1e-4:0.1;
x = (1+cos(2*pi*50*t)).*cos(2*pi*1000*t);

plot(t,x)
xlim([0 0.04]);

clf;
y = hilbert(x);
env = abs(y);
plot_param = {'Color', [0.6 0.1 0.2],'Linewidth',2}; 

plot(t,x)
hold on
plot(t,[-1]*env,plot_param{:})
plot(t,[1]*env,plot_param{:})
plot(t,[-1;1]*env,plot_param{:})
hold off
xlim([0 0.04])
title('Hilbert Envelope')

%**************** adjust Hilbert filter length
clear;

fl1 = 12;
[up1,lo1] = envelope(x,fl1,'analytic');
fl2 = 30;
[up2,lo2] = envelope(x,fl2,'analytic');
param_small = {'Color',[0.9 0.4 0.1],'Linewidth',2};
param_large = {'Color',[0 0.4 0],'Linewidth',2};

plot(t,x)
hold on
p1 = plot(t,up1,param_small{:});
plot(t,lo1,param_small{:});
p2 = plot(t,up2,param_large{:});
plot(t,lo2,param_large{:});
hold off

legend([p1 p2],'fl = 12','fl = 30')
xlim([0 0.04])
title('Analytic Envelope')

%**************** adjust length of window
wl1 = 3;
[up1,lo1] = envelope(x,wl1,'rms');
wl2 = 5;
[up2,lo2] = envelope(x,wl2,'rms');
wl3 = 300;
[up3,lo3] = envelope(x,wl3,'rms');

plot(t,x)
hold on
p1 = plot(t,up1,param_small{:});
plot(t,lo1,param_small{:});
p2 = plot(t,up2,plot_param{:});
plot(t,lo2,plot_param{:});
p3 = plot(t,up3,param_large{:});
plot(t,lo3,param_large{:})
hold off

legend([p1 p2 p3],'wl = 3','wl = 5','wl = 300')
xlim([0 0.04])
title('RMS Envelope')

%********** spline interpolation over local maxima
np1 = 5;
[up1,lo1] = envelope(x,np1,'peak');
np2 = 50;
[up2,lo2] = envelope(x,np2,'peak');

plot(t,x)
hold on
p1 = plot(t,up1,param_small{:});
plot(t,lo1,param_small{:})
p2 = plot(t,up2,param_large{:});
plot(t,lo2,param_large{:})
hold off

legend([p1 p2],'np = 5','np = 50')
xlim([0 0.04])
title('Peak Envelope')

%% https://www.mathworks.com/help/signal/ref/envelope.html
clear;
t = 0:1/2000:2-1/2000;
q = chirp(t-2,4,1/2,6,'quadratic',100,'convex').*exp(-4*(t-1).^2);
plot(t,q)

[up,lo] = envelope(q);
hold on
plot(t,up,t,lo,'linewidth',1.5)
legend('q','up','lo')
hold off
%***
clear;
t = 0:1/2000:2-1/2000;
q = chirp(t-2,4,1/2,6,'quadratic',100,'convex').*exp(-4*(t-1).^2);
envelope(q)



