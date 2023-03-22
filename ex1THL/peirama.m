%TELECOMMUNICATION SYSTEMS I
%exercise 1
%authors-Panagiotis Sklavos/Giannis Peridis

%MATLAB CODE
clear all;
close all;
clc;
format compact;

disp('press any key to see the two Rff graphs of the theoritical questions')
pause
%Th.1 Let T=10. ,Ts =0.01
T=10;
Ts=0.01;
t=-25:Ts:25;
phi=zeros(1,length(t));
phi(abs(t)<= T/2)=1/sqrt(T);

%plot(t,phi)
t_rev=-t(end:-1:1);
phi_rev=phi(end:-1:1);
tconv=t(1)+t_rev(1):Ts:t(end)+t_rev(end);
Rff=conv(phi,phi_rev)*Ts;
plot(tconv,Rff)
xlabel('\tau');
ylabel('Rff');
title('Rff of phi \Theta.1');


%Th.3 Let T=10. ,Ts =0.01
figure();
Rff=zeros(1,length(tconv));
Rff(-T<=tconv & tconv<-T/2)=-1-tconv(-T<=tconv & tconv<-T/2)/T;
Rff(-T/2<=tconv & tconv<0)=1+3*tconv(-T/2<=tconv & tconv<0)/T;
Rff(0<=tconv & tconv<T/2)=1-3*tconv(0<=tconv & tconv<T/2)/T;
Rff(T/2<=tconv & tconv<T)=-1+tconv(T/2<=tconv & tconv<T)/T;
plot(tconv,Rff)
xlabel('\tau');
ylabel('Rff');
title('Rff of phi \Theta.3');

%A.1
%intialization of symbol period T,oversampling factor over,
%parameter A(half duration of the pulse) and roll-off factor a
T=10^-3;
over=10;
A=4;
a=[0 0.5 1];

%creating srrc pulses
[f1,t1]=srrc_pulse(T,over,A,a(1));
[f2,t2]=srrc_pulse(T,over,A,a(2));
[f3,t3]=srrc_pulse(T,over,A,a(3));

%creating graph of the pulses in the same plot
disp('press any key to see the graph of the srrc pulses with diferrent a')
pause
figure();
plot(t1,f1);
hold on;
plot(t2,f2);
plot(t3,f3);
title('srrc pulses with different a')
xlabel('t=time in seconds');
ylabel('f=srrc pulse');
legend('a = 0', 'a = 0.5', 'a = 1');

%A.2
%intialization of sampling period Ts,sampling frequency Fs,
%parameter Nf(number of equidistant points) and F axis
Ts=T/over;
Fs=1/Ts;
Nf=2048;
F=[-Fs/2:Fs/Nf:Fs/2-Fs/Nf];

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
figure();
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
figure();
semilogy(F,edsF1);
hold on;
semilogy(F,edsF2);
semilogy(F,edsF3);
title('energy density spectrums with logarithmic yaxis');
xlabel('F=frequency in Hz');
ylabel('fedF=energy density spectrum of srrc pulses');
legend('a = 0', 'a = 0.5', 'a = 1');

%A.3
%creating the band windth of energy density spectrum in theory
disp(' ')
disp('the band windth of enery density spectrums in theory')
disp(' ')
BWedsF1=(1+a(1))/(2*T)
BWedsF2=(1+a(2))/(2*T)
BWedsF3=(1+a(3))/(2*T)

%creating two horizontical lines with diffrent values c1 and c2
c1=T/10^3;
for ki= 1:length(F)-1
    c1=[c1 T/10^3];
end
c2=T/10^5;
for ki= 1:length(F)-1
    c2=[c2 T/10^5];
end

disp('press any key to see the same graph ,with two constant values c1 and c2')
pause
semilogy(F,c1);
semilogy(F,c2);
legend('a = 0', 'a = 0.5', 'a = 1','c1','c2');

%B.1 
%intialization of parameter A and vector of integer values k
A = 5;
k = 0:2*A;
%creation of srrc pulses
[f1,t1]=srrc_pulse(T,over,A,a(1));
[f2,t2]=srrc_pulse(T,over,A,a(2));
[f3,t3]=srrc_pulse(T,over,A,a(3));
%initialization of integrals to vectos filled with zeros
integral1=zeros(1,length(k));
integral2=integral1;
integral3=integral1;

%for loop for every value of vector k
for ki=0:length(k)-1
    
   %creation of the delayed srrc pulse for a=0
    f1kt = [zeros(1, over*ki) f1(1:end-ki*over)];
   %creation of product (srrc pulse)*(srrc pulse delayed by kT) for a=0
    product1 = f1.*f1kt;
    if(ki<4)
        disp('press any key to see the graphs of (srrc pulse) and (srrc pulse delayed by kT)for a=0:')
        pause
        figure()
        plot(t1,f1)
        hold on
        plot(t1,f1kt);
        title('for a=0:(srrc pulse) and (srrc pulse delayed by kT)')
        ylabel('\phi(t),\phi(T-kT)')
        xlabel('t=time in seconds')
        hold off
        disp('press any key to see the graph of(srrc pulse)*(srrc pulse delayed by kT) for a=0:')
        pause
        figure()
        plot(t1,product1)
        title('for a=0:(srrc pulse)*(srrc pulse delayed by kT)')
        ylabel('\phi(t),\phi(T-kT)')
        xlabel('t=time in seconds')
        %creation of integral of (srrc pulse)*(srrc pulse delayed by kT)
        %for a=0
        integral1(ki+1)=sum(product1)*Ts;
    end
    
   %creation of the delayed srrc pulse for a=0.5
   f2kt = [zeros(1, over*ki) f2(1:end-ki*over)];
   %creation of product (srrc pulse)*(srrc pulse delayed by kT) for a=0.5
    product2 = f2.*f2kt;
    if(ki<4)
        disp('press any key to see the graphs of (srrc pulse) and (srrc pulse delayed by kT)for a=0.5:')
        pause
        figure()
        plot(t2,f2)
        hold on
        plot(t2,f2kt);
        title('for a=0.5:(srrc pulse) and (srrc pulse delayed by kT)')
        ylabel('\phi(t),\phi(T-kT)')
        xlabel('t=time in seconds')
        hold off
        disp('press any key to see the graph of(srrc pulse)*(srrc pulse delayed by kT) for a=0.5:')
        pause
        figure()
        plot(t2,product2)
        title('for a=0.5:(srrc pulse)*(srrc pulse delayed by kT)')
        ylabel('\phi(t),\phi(T-kT)')
        xlabel('t=time in seconds')
        %creation of integral of (srrc pulse)*(srrc pulse delayed by kT)
        %for a=0.5
        integral2(ki+1)=sum(product2)*Ts;
    end
    
    %creation of the delayed srrc pulse for a=1
    f3kt = [zeros(1, over*ki) f3(1:end-ki*over)];
    %creation of product (srrc pulse)*(srrc pulse delayed by kT) for a=1
    product3 = f3.*f3kt;
    if(ki<4)
        disp('press any key to see the graphs of (srrc pulse) and (srrc pulse delayed by kT)for a=1:')
        pause
        figure()
        plot(t3,f3)
        hold on
        plot(t3,f3kt);
        title('for a=1:(srrc pulse) and (srrc pulse delayed by kT)')
        ylabel('\phi(t),\phi(T-kT)')
        xlabel('t=time in seconds')
        hold off
        disp('press any key to see the graph of(srrc pulse)*(srrc pulse delayed by kT) for a=1:')
        pause
        figure()
        plot(t3,product3)
        title('for a=1:(srrc pulse)*(srrc pulse delayed by kT)')
        ylabel('\phi(t),\phi(T-kT)')
        xlabel('t=time in seconds')
        %creation of integral of (srrc pulse)*(srrc pulse delayed by kT)
        %for a=1
        integral3(ki+1)=sum(product3)*Ts;
    end
end

disp(' ')
disp('the values of the integral of (srrc pulse)*(srrc pulse delayed by kT) for a=0,0.5,1 and k=0,1,2,3 ')
disp(' ')
integral1(1:4)
integral2(1:4)
integral3(1:4)

%C.1
%initialization of symbol period T,oversampling factor over,roll-off factor
%a,parameter A(half duration of the pulse),N number of bits
T=0.1;
over=10;
a=0.5;
A=5;
N=100;
T_s=T/over;

%call of the created function that transfers bits to 2PAM
disp('a random sequence of bits')
b=(sign(randn(N,1))+1)/2;

%C.2
%a
%call of the created function that transfers bits to 2PAM
disp('the transformed 2PAM sequence')
X=bits_to_2PAM(b);

%C.2.b
%creation of x_delta
X_delta=1/T_s*upsample(X,over);
t=(0:T_s:N*T-T_s);
%creation of the graph of x_delta
disp('press any key to see the graph of x_delta')
pause
figure()
plot(t,X_delta);
title('X_delta pulse')
xlabel('X_delta(t)')
ylabel('t=time in seconds')

%c
%creating an srrc pulse
[f4,t4]=srrc_pulse(T,over,A,a);
%creating th time axis for the convolution of f4 and X_delta
tconv=[t(1)+t4(1):T_s:t(end)+t4(end)];
%convolution of f4 and X_delta
X_t=conv(f4,X_delta)*T_s;
disp('press any key to see the convolution of f4 and X_delta')
pause
figure()
plot(tconv,X_t);
title('convolution of \phi(t) and X_delta')
ylabel('convolution of \phi(t) and X_delta')
xlabel('t=time in seconds')

%d
%convolution of f4 and X_t
%no need to change f4 because it is symmetrical to the time 
%axis(\phi(t)=\phi(-t))
Z=conv(f4,X_t)*T_s;
%creating th time axis for the convolution of f4 and X_t
tconv2=[tconv(1)+t4(1):T_s:tconv(end)+t4(end)];
disp('press any key to see the convolution of \phi(-t) and X_t')
pause
figure()
hold on
plot(tconv2,Z);
stem([0:N-1]*T,X);
title('convolution of \phi(-t) and X_t')
ylabel('convolution of \phi(-t) and X_t')
xlabel('t=time in seconds')
legend('Z(t)','X')


%Printing the values of Z for kt
for i = 1:N
    values = 'Z(%d)=%8f    ||    X(%d)=%d     ||    k=%d \n';
    fprintf(values,100+i*10,Z(90+i*10),i-1,X(i),i-1)
end
