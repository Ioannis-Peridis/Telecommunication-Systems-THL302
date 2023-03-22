function [est_X,est_bit_seq] = detect_PSK_8(Y)

% This function takes the sequence that reaches the receiver, and tries to 
% determine what the initial sequence of symbols was, and subsequently,
% what was the initial bit sequence.
% INPUT:
% Y: The sequence at the receiver
% OUTPUT:
% est_X: The estimated sequence of symbols the transmitter sent
% est_bit_seq: The estimated sequence of bits the transmitter sent

% Initializations
len_Y = length(Y);
est_X = zeros(len_Y,2);
est_bit_seq = zeros(3*len_Y,1);
len_Y = length(Y);


% All the possible positions on the complex circle:
zero = [cos(2*pi*0/8), sin(2*pi*0/8)];
one = [cos(2*pi*1/8), sin(2*pi*1/8)];
two = [cos(2*pi*2/8), sin(2*pi*2/8)];
three = [cos(2*pi*3/8), sin(2*pi*3/8)];
four = [cos(2*pi*4/8), sin(2*pi*4/8)];
five = [cos(2*pi*5/8), sin(2*pi*5/8)];
six = [cos(2*pi*6/8), sin(2*pi*6/8)];
seven = [cos(2*pi*7/8), sin(2*pi*7/8)];

positions = [zero; one; two; three; four; five; six; seven];
len_pos = length(positions);

for i=1: len_Y
    best_est = zero;
    min_dist = 100000;
    best_est_index = 0;
    
    for j=1: len_pos
        if(norm(Y(i,:)-positions(j,:),2) < min_dist)
            min_dist = norm(Y(i,:)-positions(j,:),2);
            best_est = positions(j,:);
            best_est_index = j;
        end
    end
    est_X(i,1) = best_est(1);
    est_X(i,2) = best_est(2);
    
    % Now the bit sequence for the specific symbol is generated
    % position:  0  ->  1  ->  2  ->  3  ->  4  ->  5  ->  6  ->  7
    % in index:  1  ->  2  ->  3  ->  4  ->  5  ->  6  ->  7  ->  8
    % bits    : 000 -> 001 -> 011 -> 010 -> 110 -> 111 -> 101 -> 100
  
    % The index in the symbol sequence.
    index = 3*(i-1)+1;
    if(best_est_index == 1)
       est_bit_seq(index) = 0;
       est_bit_seq(index+1) = 0;
       est_bit_seq(index+2) = 0;
    elseif(best_est_index == 2)
       est_bit_seq(index) = 0;
       est_bit_seq(index+1) = 0;
       est_bit_seq(index+2) = 1;
    elseif(best_est_index == 3)
       est_bit_seq(index) = 0;
       est_bit_seq(index+1) = 1;
       est_bit_seq(index+2) = 1;
    elseif(best_est_index == 4)
       est_bit_seq(index) = 0;
       est_bit_seq(index+1) = 1;
       est_bit_seq(index+2) = 0;
    elseif(best_est_index == 5)
       est_bit_seq(index) = 1;
       est_bit_seq(index+1) = 1;
       est_bit_seq(index+2) = 0;
    elseif(best_est_index == 6)
       est_bit_seq(index) = 1;
       est_bit_seq(index+1) = 1;
       est_bit_seq(index+2) = 1;
    elseif(best_est_index == 7)
       est_bit_seq(index) = 1;
       est_bit_seq(index+1) = 0;
       est_bit_seq(index+2) = 1;
    elseif(best_est_index == 8)
       est_bit_seq(index) = 1;
       est_bit_seq(index+1) = 0;
       est_bit_seq(index+2) = 0;       
    end
end
end

