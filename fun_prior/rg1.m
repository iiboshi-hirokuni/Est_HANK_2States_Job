function x = rg1(a)

if a > 1;
  x = rg2(a);
elseif a < 1;
  a = a + 1;
  u = rand();
  x = rg2(a)*u^(1/a);
elseif a == 1;
  x = -log(rand());
end;
