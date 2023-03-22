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