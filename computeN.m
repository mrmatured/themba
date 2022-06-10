function N =  computeN(psi1, psi2, L, N)
  N0 = trapz(abs(psi1))*L/N
  N1 = trapz(abs(psi2))*L/N
  N  = abs((N1 - N0)/N0);
end

