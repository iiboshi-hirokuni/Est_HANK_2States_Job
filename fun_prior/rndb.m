function x = rndb(a,b)
% beta

a1n = rg1(a);
a1d = rg1(b);
x = a1n / (a1n + a1d);
