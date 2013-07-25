function nextstate = RK4Integrator(DerivFcn,t,state,dt)
%
% RK4Integrator
%
% DerivFcn  - Fcn Handle to eqns of motion: d/dt(state) = DerivFcn(t,state)
% t         - current time
% state     - current state
% dt        - time interval
% nextstate - state at time t+dt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth Order Runge Kutta Integrator
%
% This integration routine propagates the state forward in time by using
% four approximations of the state derivative and then combining them.
%
% This integrator is kinda awesome. You can use larg(er) timesteps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


statedot1 = feval(DerivFcn,t,state);
statedot2 = feval(DerivFcn,t+dt/2,state+statedot1*dt/2);
statedot3 = feval(DerivFcn,t+dt/2,state+statedot2*dt/2);
statedot4 = feval(DerivFcn,t+dt,state+statedot3*dt);

statedot = 1/6*(statedot1 + 2*statedot2 + 2*statedot3 + statedot4);

nextstate = state + statedot*dt;