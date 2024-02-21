function [vect] = fun_std2comp(tens)
%fun_std2comp translates 2nd order tensor given in standard notation (symm
%3x3 matrix) to compressed notation (1x6 vector)
%   Input:  tens ... 2nd order tensor as 3x3 matrix
%   Output: vect ... 2nd order tensor as 1x6 matrix


% second order unity tensor
vect=[tens(1,1),tens(2,2),tens(3,3),sqrt(2)*tens(3,2),sqrt(2)*tens(1,3),sqrt(2)*tens(2,1)];
end

