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
