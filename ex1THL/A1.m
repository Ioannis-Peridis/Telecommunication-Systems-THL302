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
figure();
plot(t1,f1);
hold on;
plot(t2,f2);
plot(t3,f3);
title('srrc pulses with different a')
xlabel('t=time in seconds');
ylabel('f=srrc pulse');
legend('a = 0', 'a = 0.5', 'a = 1');

