function [tens] = fun_comp2std(vect)
%fun_comp2std translates 2nd order tensor given in compressed notation (1x6
%vector) to standard notation (symm 3x3 matrix) 
%   Input:  vect ... 2nd order tensor as 1x6 matrix
%   Output: tens ... 2nd order tensor as 3x3 matrix

% second order unity tensor
tens=[vect(1),vect(6)/sqrt(2),vect(5)/sqrt(2);...
      vect(6)/sqrt(2),vect(2),vect(4)/sqrt(2);...
      vect(5)/sqrt(2),vect(4)/sqrt(2),vect(3)];
end

