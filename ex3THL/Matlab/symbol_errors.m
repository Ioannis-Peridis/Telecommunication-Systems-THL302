function [num_of_symbol_errors] = symbol_errors(est_X,X)

% This function takes as input the sequence that reaches the receiver, 
% and the one that the transmitter sent, and attempts to find the number 
% of mistakes that have occured in the process.
% INPUT:
% est_X: The estimated sequence of symbols, the transmitter propably sent
% X    : The real sequence of symbols, sent by the transmitter.
% OUTPUT:
% num_of_symbol_errors: The number of errors that happened.

% Initializations
num_of_symbol_errors = 0;
len_X = length(X);

for i=1: len_X
  if(est_X(i) ~= X(i))
    num_of_symbol_errors = num_of_symbol_errors + 1;
  end
end

