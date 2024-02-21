function [I2] = fun_I2()
%fun_I2 provides access to second order unity tensors
%   no input, output is a structure

% second order unity tensor
I2.I=[1,0,0;0,1,0;0,0,1];
I2.mandel=[1;1;1;0;0;0];
end

