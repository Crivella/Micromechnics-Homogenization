function [Ttrans] = fun_trans(T,azi,zeni)
Q4=fun_Q4_bp(azi,zeni);
Q4t=transpose(Q4);
Ttrans=Q4*T*Q4t;
end

