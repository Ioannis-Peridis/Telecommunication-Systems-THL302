function [X] = bits_to_PSK_8(bit_seq)

% This function takes as input a sequence of bits, and converts it into
% a sequence of 8-PSK symbols by using Gray coding. Specifically the 
% coding is: 000 -> 001 -> 011 -> 010 -> 110 -> 111 -> 101 -> 100
% position :  0  ->  1  ->  2  ->  3  ->  4  ->  5  ->  6  ->  7
% Input: The bit sequence that is coded.
% Output: The 8-PSK sequence

%Initializations
N = length(bit_seq)/3;
X = zeros(N,2);

for k = 1: 3: size(bit_seq)
    % The index in the symbol sequence.
    index = (k-1)/3 +1;
    
    % Check bits by 3, and create symbols as the the Gray code implies.
    if(bit_seq(k) == 0 && bit_seq(k+1) == 0 && bit_seq(k+2) == 0)
        pos = 0;
    elseif(bit_seq(k) == 0 && bit_seq(k+1) == 0 && bit_seq(k+2) == 1)
        pos = 1;
    elseif(bit_seq(k) == 0 && bit_seq(k+1) == 1 && bit_seq(k+2) == 1)
        pos = 2;
    elseif(bit_seq(k) == 0 && bit_seq(k+1) == 1 && bit_seq(k+2) == 0)
        pos = 3;
    elseif(bit_seq(k) == 1 && bit_seq(k+1) == 1 && bit_seq(k+2) == 0)
        pos = 4;
    elseif(bit_seq(k) == 1 && bit_seq(k+1) == 1 && bit_seq(k+2) == 1)
        pos = 5;
    elseif(bit_seq(k) == 1 && bit_seq(k+1) == 0 && bit_seq(k+2) == 1)
        pos = 6;
    elseif(bit_seq(k) == 1 && bit_seq(k+1) == 0 && bit_seq(k+2) == 0)
        pos = 7;
    end
    X(index,1) = cos(2*pi*pos/8);
    X(index,2) = sin(2*pi*pos/8);
end   
end