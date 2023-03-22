%A.2
%intialization of sampling period Ts,sampling frequency Fs,
%parameter Nf(number of equidistant points) and F axis
Ts=T/over;
Fs=1/Ts;
Nf=2048;
F=-Fs/2:Fs/Nf:Fs/2-Fs/Nf;
%creating fast fourier transformations
F1= fftshift(fft(f1,Nf)*Ts);
F2= fftshift(fft(f2,Nf)*Ts);
F3= fftshift(fft(f3,Nf)*Ts);
%creating the energy density spectrumd
edsF1=abs(F1).^2;
edsF2=abs(F2).^2;
edsF3=abs(F3).^2;
%creating graph of the energy density spectrums
%in the same plot
disp('press any key to see the graph of the energy density spectrums')
pause
figure(2);
plot(F,edsF1);
hold on;
plot(F,edsF2);
plot(F,edsF3);
title('energy density spectrums');
xlabel('F=frequency in Hz');
ylabel('fedF=energy density spectrum of srrc pulses');
legend('a = 0', 'a = 0.5', 'a = 1');
%creating graph of the energy density spectrums
%in the same semilogy
disp('press any key to see the logarithmic yaxis graph of the energy density spectrums')
pause
figure(3);
semilogy(F,edsF1);
hold on;
semilogy(F,edsF2);
semilogy(F,edsF3);
title('energy density spectrums with logarithmic yaxis');
xlabel('F=frequency in Hz');
ylabel('fedF=energy density spectrum of srrc pulses');
legend('a = 0', 'a = 0.5', 'a = 1');