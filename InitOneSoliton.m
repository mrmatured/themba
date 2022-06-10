function [u,v] = InitOneSoliton(x, off, v1)
  u = sech(x-off).*exp(1i*v1*(x-off));
  v = sech(x-off).*exp(1i*v1*(x-off));
end

