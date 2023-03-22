%TELECOMMUNICATION SYSTEMS I
%exercise 3
%authors-Panagiotis Sklavos/Giannis Peridis

%MATLAB CODE
clear all;
close all;
clc;
format compact;

%=====================================================================================
%A.1

%N number of bits
N=100;
%creating a'a random sequence of 3*N independent and  equally possible bits
%disp('a random sequence of 3*N independent and  equally possible bits')
b=(sign(randn(3*N,1))+1)/2;

%=====================================================================================
%A.2

%call of the created function that transfers bits to 8PSK
X = bits_to_PSK_8(b);
Xi = transpose(X(:,1));
Xq = transpose(X(:,2));
%creation of the graph of XQ and XI symbols sequencies in discrete time in
figure();
subplot(2, 1, 1)
stem(Xi)
title('X_I symbols sequency');
xlabel('DISCRETE TIME');
ylabel('X_I');
subplot(2, 1, 2)
stem(Xq)
title('X_Q symbols sequency');
xlabel('DISCRETE TIME');
ylabel('X_Q');

%=====================================================================================
%A.3

%intialization of symbol period T,oversampling factor over,
%parameter A(half duration of the pulse)  and roll-off factor  
T=10^-3;
over=10;
A=4;
a=0.5;
%creating srrc pulses
[phi,t]=srrc_pulse(T,over,A,a);

%intialization of sampling period Ts,sampling frequency Fs,
%parameter Nf(number of equidistant points) and F axis
Ts=T/over;
Fs=1/Ts;
Nf=2048;
F=[-Fs/2:Fs/Nf:Fs/2-Fs/Nf];

%time axis of Xi and Xq
t1=(0:Ts:N*T-Ts);

%Creation of continious Xi and Xq:

%creation of X_delta_i and X_delta_q
X_delta_i=1/Ts*upsample(Xi,over);
X_delta_q=1/Ts*upsample(Xq,over);
%creation of Xic=sum(deltai*phi(t-n*T)) and Xiq=sum(deltaq*phi(t-n*T))
%convolution of phi and X_delta_i and X_delta_q
tconv=[t(1)+t1(1):Ts:t(end)+t1(end)];
Xic=conv(phi,X_delta_i)*Ts;
Xqc=conv(phi,X_delta_q)*Ts;

%creation of the graph of XIc and XQc symbols sequencies in continious time in
figure();
subplot(2,1,1)
plot(tconv,Xic)
title('X_Ic symbols sequency');
xlabel('CONTINIOUS TIME');
ylabel('X_I');
subplot(2,1,2)
plot(tconv,Xqc)
title('X_Qc symbols sequency');
xlabel('CONTINIOUS TIME');
ylabel('X_Q');

%creating fast fourier transformation of Xic and Xqc
XIC=fftshift(fft(Xic,Nf)*Ts);
XQC=fftshift(fft(Xqc,Nf)*Ts);
%the total duration time
Ttotal=length(tconv)*Ts;

%creation of a periodogram of one of the implementations of XIC and XIQ
P_xic=abs(XIC).^2/Ttotal;
P_xqc=abs(XQC).^2/Ttotal;

%creation of periodogram graph of XQc and XIc in logarithmic y axis 
figure();
subplot(2,1,1)
semilogy(F,P_xic)
title('Periodogram of X_Ic ');
xlabel('F=frequency in Hz');
ylabel('P(X_Ic)');
subplot(2,1,2)
semilogy(F,P_xqc)
title('Periodogram of X_Qc ');
xlabel('F=frequency in Hz');
ylabel('P(X_Qc)');

%=====================================================================================
%A4

%initialization of carrier frequency Fo
Fo=2000;
%creation of the carrier
z1=2*cos(2*pi*Fo*tconv);
z2=-2*sin(2*pi*Fo*tconv);
%creation of modulated signals XIt and XQt
Xit=Xic.*z1;
Xqt=Xqc.*(z2);

%creation of the graphs of the modulated signals XIt and XQt 
figure();
subplot(2,1,1)
plot(tconv, Xit)
title('XIc modulated by carrier cos');
xlabel('t=time in seconds');
ylabel('X_It');
subplot(2,1,2)
plot(tconv, Xqt)
title('XQc modulated by carrier -sin');
xlabel('t=time in seconds');
ylabel('X_Qt');

%creating fast fourier transformation of Xit and Xqt
XIT=fftshift(fft(Xit,Nf)*Ts);
XQT=fftshift(fft(Xqt,Nf)*Ts);
%the total duration time
Ttotal=length(tconv)*Ts;

%creation of a periodogram of one of the implementations of XIt and XQt
P_xit=abs(XIT).^2/Ttotal;
P_xqt=abs(XQT).^2/Ttotal;

%creation of periodogram graph of XQt and XIt in logarithmic y axis 
figure();
subplot(2,1,1)
semilogy(F,P_xit)
title('periodogramm of X_It ');
xlabel('F=frequency in HZ');
ylabel('P(X_It)');
subplot(2,1,2)
semilogy(F,P_xqt)
title('periodogramm of X_Qt ');
xlabel('F=frequency in HZ');
ylabel('P(X_Qt)');

%=====================================================================================
%A.5,6

%creation of the channel's input X(t)
X_t=Xit+Xqt;

%creation of the inputs graph 
figure();
plot(tconv,X_t)
title('channels input WITHOUT NOISE FACTOR ');
xlabel('t=time in seconds');
ylabel('X(t)');

%creating fast fourier transformation of Xic and Xqc
XF=fftshift(fft(X_t,Nf)*Ts);
%the total duration time
Ttotal=length(tconv)*Ts;
%creation of a periodogram of one of the implementations of X(t)
P_x=abs(XF).^2/Ttotal;

%creation of periodogram graph of X(t) in logarithmic y axis
figure();
semilogy(F,P_x)
title('periodogramm of X(t) ');
xlabel('F=frequency in HZ');
ylabel('P(X)');

%=====================================================================================
%A.7

%insertion of white gaussian noise
SNR_DB=20;
var = sqrt(1/(Ts*10^(SNR_DB/10)));
W_t = var*randn(length(X_t),1);
Y_t = zeros(1, length(X_t));

for i=1 : length(X_t)
    Y_t(i) = X_t(i) + W_t(i);
end
%creation of X_t with noise factor
figure();
plot(tconv,Y_t)
title('channels input WITH NOISE FACTOR ');
xlabel('t=time in seconds');
ylabel('X(t) with noise');

%creation of the difference of X_t with and without noise 
diff=Y_t-X_t;
%creation of the difference graph
figure();
plot(tconv,diff)
title('difference of the input signals with and without noise');
xlabel('t=time in seconds');
ylabel('X(t)with noise-X(t)without noise ');

%=====================================================================================
%A.8

%creation of remodulated signals Yz1 and Yz2
Yz1=Y_t.*cos(2*pi*Fo*tconv);
Yz2=Y_t.*(-sin(2*pi*Fo*tconv));

%creation of ther graph of Yz1 and Yz2
figure();
subplot(2,1,1)
plot(tconv,Yz1)
title('Y_t remodulated by carrier cos ');
xlabel('t=time in seconds');
ylabel('Yz1');
subplot(2,1,2)
plot(tconv,Yz2)
title('Y_t remodulated by carrier -sin');
xlabel('t=time in seconds');
ylabel('Yz2');

%creating fast fourier transformation of Yz1 and Yz2   
YZ1=fftshift(fft(Yz1,Nf)*Ts);
YZ2=fftshift(fft(Yz2,Nf)*Ts);
%the total duration time
Ttotal=length(tconv)*Ts;
%creation of a periodogram of one of the implementations of Yz1 and Yz2
P_yz1=abs(YZ1).^2/Ttotal;
P_yz2=abs(YZ2).^2/Ttotal;

%creation of the graph of periodogramms of YZ1 and YZ2
figure();
subplot(2,1,1)
semilogy(F,P_yz1)
title('periodogramm of the  Y_t remodulated by carrier cos ');
xlabel('F=frequency in HZ');
ylabel('P(YZ1)');
subplot(2,1,2)
semilogy(F,P_yz2)
title('periodogramm of the Y_t remodulated by carrier -sin ');
xlabel('F=frequency in HZ');
ylabel('P(YZ2)');

figure();
subplot(2,1,1)
plot(F,P_yz1)
title('periodogramm of the  Y_t remodulated by carrier cos ');
xlabel('F=frequency in HZ');
ylabel('P(YZ1)');
subplot(2,1,2)
plot(F,P_yz2)
title('periodogramm of the Y_t remodulated by carrier -sin ');
xlabel('F=frequency in HZ');
ylabel('P(YZ2)');


%=====================================================================================
%A9

%creation of the new time axis 
new_tconv=tconv(1)+t(1):Ts:tconv(end)+t(end);

%Yz1new=the new Xit after the filter
Yz1new=conv(Yz1,phi)*Ts;
%Yz2new=the new Xqt after the filter
Yz2new=conv(Yz2,phi)*Ts;

%creation of the new waveforms Yz1new and Yz2new
figure();
subplot(2,1,1)
plot(new_tconv,Yz1new)
title('Yz1 after the filter= Xit*');
xlabel('t=time in seconds');
ylabel('Xit*');
subplot(2,1,2)
plot(new_tconv,Yz2new)
title('Yz2 after the filter= Xqt*');
xlabel('t=time in seconds');
ylabel('Xqt*');

%creating fast fourier transformation of Yz1new and Yz2new   
YZ1NEW=fftshift(fft(Yz1new,Nf)*Ts);
YZ2NEW=fftshift(fft(Yz2new,Nf)*Ts);
%the total duration time
Ttotal=length(tconv)*Ts;
%creation of a periodogram of one of the implementations of Yz1 and Yz2
P_yz1new=abs(YZ1NEW).^2/Ttotal;
P_yz2new=abs(YZ2NEW).^2/Ttotal;

%creation of the graph of periodogramms of YZ1NEW and YZ2NEW
figure();
subplot(2,1,1)
semilogy(F,P_yz1new)
title('periodogramm of Xit* ');
xlabel('F=frequency in HZ');
ylabel('P(XIF*)');
subplot(2,1,2)
semilogy(F,P_yz2new)
title('periodogramm of Xqt* ');
xlabel('F=frequency in HZ');
ylabel('P(XQF*)');

figure();
subplot(2,1,1)
plot(F,P_yz1new)
title('periodogramm of Xit* with plot');
xlabel('F=frequency in HZ');
ylabel('P(XIF*)');
subplot(2,1,2)
plot(F,P_yz2new)
title('periodogramm of Xqt* ');
xlabel('F=frequency in HZ');
ylabel('P(XQF*)');
%=====================================================================================    
%A10

%sampling of the output of the filters , Xit* and Xqt*
n=0:N-1; 
yi(n+1) = Yz1new(n*over + 80); 
yq(n+1) = Yz2new(n*over + 80);
Y = [transpose(yi),transpose(yq)];

%creation of the new sampled signals with scatter graph     
scatterplot(Y)
title("the received sampled symbols")

scatterplot(X)
title("the initial symbols")
%=====================================================================================
%A11

Y=[transpose(yi),transpose(yq)];
%calling of the detect function
[est_X, est_bit_seq] = detect_PSK_8(Y);

%=====================================================================================
%A12

%calling of the symbol errors function
num_of_symbol_errors = symbol_errors(est_X,X);
fprintf("number of error symbols is %d \n", num_of_symbol_errors);

%=====================================================================================   
%A13

%calling of the bit errors function
num_of_bit_errors = bit_errors(est_bit_seq,b);
fprintf("number of error bits is %d", num_of_bit_errors);

%=====================================================================================    
%B

%number of times the process is runed
K=1000;
%SNR in decibels
SNR_db= [-2:2:16];
%the monte carlo probabilities
%initializing the symbol error probability
P_Esymbol=zeros(1,length(SNR_db));
%initializing the bit error probability
P_Ebit=zeros(1,length(SNR_db));
%initializing the upper and lower bounds
upper_bound = zeros(1,length(SNR_db));
lower_bound = zeros(1,length(SNR_db));

for i=1:length(SNR_db)
    
    %initializing to zero the number of symbol and bit errors
    num_of_symbol_errors=0;
    num_of_bit_errors=0;
    
    for j=1:K
        
        %creation of bit sequence
        b=(sign(randn(3*N,1))+1)/2;
        
        %call of the created function that transfers bits to 8PSK
        X = bits_to_PSK_8(b);
        Xi = transpose(X(:,1));
        Xq = transpose(X(:,2));
        
        %creating the srrc pulse
        [phi,t]=srrc_pulse(T,over,A,a);
        
        %time axis of Xi and Xq
        t1=(0:Ts:N*T-Ts);
        %creation of X_delta_i and X_delta_q
        X_delta_i=1/Ts*upsample(Xi,over);
        X_delta_q=1/Ts*upsample(Xq,over);
        
        %convolution of phi and X_delta_i and X_delta_q
        tconv=[t(1)+t1(1):Ts:t(end)+t1(end)];
        Xic=conv(phi,X_delta_i)*Ts;
        Xqc=conv(phi,X_delta_q)*Ts;
        
        %creation of the carrier
        z1=2*cos(2*pi*Fo*tconv);
        z2=-2*sin(2*pi*Fo*tconv);
        %creation of modulated signals XIt and XQt
        Xit=Xic.*z1;
        Xqt=Xqc.*(z2);
        
        %creation of the channel's input X(t)
        X_t=Xit+Xqt;
        
        %insertion of white gaussian noise
        disper = sqrt(1/(Ts*10^(SNR_db(i)/10)));
        W_t = disper*randn(length(X_t),1);
        Y_t = zeros(1, length(X_t));
        
        for p=1 : length(X_t)
         Y_t(p) = X_t(p) + W_t(p);
        end
        
        %creation of remodulated signals Yz1 and Yz2
        Yz1=Y_t.*cos(2*pi*Fo*tconv);
        Yz2=Y_t.*(-sin(2*pi*Fo*tconv));
        
        %creation of the new time axis 
        new_tconv=tconv(1)+t(1):Ts:tconv(end)+t(end);
        %Yz1new=the new Xit after the filter
        Yz1new=conv(Yz1,phi)*Ts;
        %Yz2new=the new Xqt after the filter
        Yz2new=conv(Yz2,phi)*Ts;

        %sampling of the output of the filters , Xit* and Xqt*
        n=0:N-1; 
        yi(n+1)=Yz1new(n*over + 80); 
        yq(n+1)=Yz2new(n*over + 80);

        %calling of the detect function
        Y=[transpose(yi),transpose(yq)];
        [est_X, est_bit_seq] = detect_PSK_8(Y);
     
        %calling of the symbol errors function and counting the total
        %number of symbol errors
        num_of_symbol_errors = num_of_symbol_errors + symbol_errors(est_X,X);

        %calling of the bit errors function and counting the total number
        %of bit errors
        num_of_bit_errors = num_of_bit_errors + bit_errors(est_bit_seq,b);
    end
    
    %calculating the probabilities of symbol and bit errors with the monte
    %carlo method
    P_Esymbol(1, i) = num_of_symbol_errors/((length(b)/3)*K);
    P_Ebit(1, i) = num_of_bit_errors/(length(b)*K);
    %calcuting the the smart upper bound for every SNR, with the given Q function
    upper_bound(i) = 2*Q(sqrt(2*(10^(SNR_db(1,i)/10)))*sin(pi/8));
    %calculating the lower bound , for every wrong symbol there is at least
    %one wrong bit
    lower_bound(i) = upper_bound(i)/3; 
 
end

%creation of the graph  of the smart upper bound and the probability of the
%symbol errors in the same plot with logarithmic y axis
figure();
semilogy(SNR_db,P_Esymbol);
hold on;
semilogy(SNR_db, upper_bound);
legend('P(Esymbols)','Smart Upper Bound');
xlabel('SNR_{db}');
ylabel('Probability of symbol error');
title('Graph with logarithmic y axis with the Theoritical and Expeimental Probability Of Symbol Error');
hold off;

%creation of the graph  of the lower bound and the probability of the
%biterrors in the same plot with logarithmic y axis
figure();
semilogy(SNR_db,P_Ebit);
hold on;
semilogy(SNR_db, lower_bound);
legend('P(Ebits)','Lower Bound');
xlabel('SNR_db');
ylabel('Probability of bit error');
title('Graph with logarithmic y axis with the Theoritical and Expeimental Probability Of Bit Error');
hold off;