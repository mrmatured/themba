%% Solving the Coupled Nonlinear Schrodinger Equation found in the context of
%% optical pulses using the Finite-Difference Method (CN) and Split-Step Method.
function [x,tdata, udata, vdata, RunTime] = Split_Step(tau,N,L,Tmax)
  global alpha x1 x2 v1 v2 off
  tic;
  h  = L/N;               % time step
  x  = [-L/2:h:L/2-h]';   % x-vector
  Nt = round(Tmax/tau);   % number of time step

  %%wave number
  k  = (2*pi/L)*[0:N/2-1 -N/2:-1]';
  k2 = k.*k;
  %% initial configuration
%%  [u,v] = InitOneSoliton(x,off, v1);          %% One Soliton
  [u, v] = Initial_Cond(x,alpha,v1,v2,x1,x2); %% two solitons
  %% computing the exponent
  EXP = exp(-1i*k2*tau);
  %%main calculation
  udata(:,1) = u; vdata(:,1) = v; t =0;tdata(1) = t;
  for tt = 2:Nt
    %update time
    t = t + tau;
    %% solution of the nonlinear part
    F = abs(u).*abs(u)-alpha*cos(abs(v).*abs(v));
    G = abs(v).*abs(v)-alpha*cos(abs(u).*abs(u));
    u = u.*exp(1i*tau*F); v = v.*exp(1i*tau*G);
    %% solution of the linear part
    u = ifft(EXP.*fft(u)); v = ifft(EXP.*fft(v));
    udata(:,tt)= u; vdata(:,tt)= v; tdata(tt)= t;
  end
  RunTime = toc
end


