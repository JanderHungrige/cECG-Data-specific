close all

   %% Time specifications:
   FS = 1000;                   % samples per second
   dt = 1/FS;                   % seconds per sample
   StopTime = 10;             % seconds
   t = (0:dt:StopTime-dt)';     % seconds
   %% Sine wave:
   Fc = 3;                     % hertz
   x = sin(2*pi*Fc*t);
   % Plot the signal versus time:
   figure;
   plot(t,x);
   ylim([-1.5 1.5]);
   xlim([0 1])
   xlabel('time');
%    title('TIme domain');
   
   x2=sin(2*pi*50*t)/5+sin(2*pi*Fc*t);
   figure;
   plot(t,x2);
   ylim([-1.5 1.5]);
   xlabel('time');
   xlim([0 1])
%    zoom xon;

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(FS*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:FS/length(x):FS/2;

figure;h2=plot(freq,10*log10(psdx));
xlabel('Frequency')
ylabel('Power/Frequency (dB/Hz)')
xlim([0 10])

N = length(x2);
xdft = fft(x2);
xdft = xdft(1:N/2+1);
psdx = (1/(FS*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:FS/length(x2):FS/2;


figure;h2=plot(freq,10*log10(psdx));
xlabel('Frequency')
ylabel('Power/Frequency (dB/Hz)')
xlim([0 60])

