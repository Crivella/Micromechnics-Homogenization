function [Trot] = fun_rot(T,azi,zeni)
Q4=fun_Q4_bp(azi,zeni);
Q4t=transpose(Q4);
Trot=Q4t*T*Q4;
end

