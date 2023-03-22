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
b=(sign(randn(N,1))+1)/2