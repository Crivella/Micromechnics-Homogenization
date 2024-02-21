function [b] = binomialkoeffizient(n,k)
if 2*k > n 
    k = n-k;
end
b = 1;
for i = 1:k
b = b * (n - k + i) / i;
end

