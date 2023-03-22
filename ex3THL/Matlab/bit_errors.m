function [num_of_bit_errors] = bit_errors(est_bit_seq,bit_seq)

% ?his function takes as input the sequence of bits that reaches the receiver, 
% and the one that the transmitter sent, and attempts to find the number 
% of errors that have occured in the process.
% INPUT:
% est_bit_seq: The estimated sequence of bits, the transmitter propably sent
% bit_seq    : The real sequence of bits, sent by the transmitter.
% OUTPUT:
% num_of_bit_errors: The number of errors that happened.

% Initializations
num_of_bit_errors = 0;
len_bit_seq = length(bit_seq);

for i=1:len_bit_seq
  if(est_bit_seq(i) ~= bit_seq(i))
    num_of_bit_errors = num_of_bit_errors + 1;
  end
end

