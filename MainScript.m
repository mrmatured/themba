%MainScript.m
% simulation the coupled system:
%                   i*u_t+u_xx +(abs(u)^2-alpha*cos(abs(v)^2))*u=0
%                   i*v_t+v_xx +(abs(v)^2-alpha*cos(abs(u)^2))*v=0
% with initial configuaration:
%           u(x,0) = sech(x-x1)*exp(1i*v1*(x-x1))+sech(x-x1)*exp(1i*v1*(x-x1))
%           v(x,0) = sech(x-x1)*exp(1i*v1*(x-x1))+sech(x-x1)*exp(1i*v1*(x-x1))
%for different values of alpha and v1, v2
%%clear all; close all; clc;
%======================== Global parameters ==================================
global alpha x1 x2 v1 v2 off
%======================== PARAMETERS =========================================
N = 1024; L =64; Tmax = 1;  tau = 0.001; alpha = 4;
%========================= Simulation Set-Up =================================
x1 = -10; x2 = 10;       % inter-soliton distance
v1 = -2; v2 = 2;         % velocity of the wave
%%v1 = -0.5; v2 = 0.5;  % velocity of the wave


%the solution via Split Step scheme
[x,tdata,  udata, vdata, RunTime] = Split_Step(tau,N,L,Tmax);

%the solution via Crank-Nicolson scheme
%%[x,tdata1, udata1, vdata1, RunTime1] = Crank_Nicolson(tau,N,L,Tmax);


%% computing the conserved quantity N via composite trapezium method
%%psi1 = abs(udata(:,1)).^2+abs(vdata(:,1)).^2; N0 = trapz(x,psi1); Nerr =[];
%%for nn = 1:16
%%  psi2 = abs(udata(:,nn)).^2+abs(vdata(:,nn)).^2; N1 = trapz(x, psi2);
%%  Nerr(nn) = abs((N1-N0)/N0);
%%end
%%max(Nerr)

%% plots of the results
[X,T]= meshgrid(x,tdata); surface(X,T,abs(udata')); colormap('jet');
colorbar("northoutside"); view(-10,75); zlabel('|u(x,t)|','fontsize',16);
xlim([-L/2 L/2]); xlabel('space:x','fontsize',16);
ylabel('time:t','fontsize',16);


%% the final profile and initial profile
figure(2)
plot(x,abs(udata(:,1)),'-b','linewidth',2,x,abs(udata(:,end)),'-k','linewidth',2)
xlabel('x','fontsize',16); ylabel('|u|','fontsize',16); grid on;
xlim([-L/2 L/2]); legend('|u(x,0)|','|u(x,t_{end})|', 'fontsize',16)
