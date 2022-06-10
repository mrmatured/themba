%% Solving the Coupled Nonlinear Schrodinger Equation found in the context of
%% optical pulses using the Finite-Difference Method (CN) and Split-Step Method.

%% Crank-Nicolson scheme for the coupled Nonlinear Schrodinger equation
function [x, tdata, udata,vdata, RunTime] = Crank_Nicolson(tau,N,L,Tmax)
  global alpha x1 x2 v1 v2 off
  tic;
  h  = L/N;             % time step
  x  = [-L/2:h:L/2-h]'; % x-vector
  Nt = round(Tmax/tau);  % number of time step
  s  = 1i*tau/(2*h^2);     % courant coefficient

  %% The nonlinear function
  F = @(uu,vv) abs(uu).*abs(uu)-alpha*cos(abs(vv).*abs(vv));
  G = @(uu,vv) abs(vv).*abs(vv)-alpha*cos(abs(uu).*abs(uu));
  %% initial configuration
%%  [u,v] = InitOneSoliton(x,off, v1);          %% One Soliton
  [u, v] = Initial_Cond(x,alpha,v1,v2,x1,x2); %% two solitons
  %% creating the A, B matrices
  LL = diag(-2*ones(N,1),0) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
  LL(1,end)=1; LL(end, 1)=1; I = eye(size(LL)); A = I-s*LL; B = I+s*LL;
  Ainv = A\eye(size(A));
  %%main calculation
  udata(:,1) = u; vdata(:,1) = v; t =0;tdata(1) = t;
  for n = 2:Nt
    %update time
    t = t + tau;
    %update solution
    F1 = F(u,v); G1 = G(u,v);
    unew = Ainv*(B*u + 1i*tau*F1.*u); vnew = Ainv*(B*v + 1i*tau*G1.*v);
    u = unew; v = vnew;
    %% solution
    udata(:,n) = u; vdata(:,n) = v;tdata(n)=t;
  end
  RunTime = toc
end


