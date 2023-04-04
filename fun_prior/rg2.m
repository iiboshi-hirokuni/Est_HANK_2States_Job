function x = rg2(a)

b = a-1;
c = 3*a - .75;
accept = 0;
while accept == 0
  u = rand();
  v = rand();
  w = u*(1-u);
  y = sqrt(c/w)*(u-.5);
  x = b+y;
  if x >= 0
    z = 64*(w^3)*(v^2);
    accept = (z <= ( 1-(2*y^2)/x ));
    if accept == 0
      accept = (log(z) <= 2*(b*log(x/b) - y));
    end
  end
end
