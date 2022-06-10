function [u, v] = Initial_Cond(x,alpha,v1,v2,x1,x2)
  u = sech(x+x1).*exp(1i*v1*(x-x1)) + sech(x-x2).*exp(1i*v2*(x-x2));
  v = sech(x+x1).*exp(1i*v1*(x-x1)) + sech(x-x2).*exp(1i*v2*(x-x2));
end

