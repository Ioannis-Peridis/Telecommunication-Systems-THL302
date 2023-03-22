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
