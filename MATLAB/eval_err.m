function [error,x] = eval_err(m,trans)
% This function simulates a first order transient based on the inertia
% value represented by 'm'. Then the simulated transient is compared to the
% real transient represented by the infromation in the struct, 'trans'. 
error = 0;

x_sim   = NaN(1,length(trans.t));
tau_sim = NaN(1,length(trans.t));

x_sim(1) = trans.x(1);

for i = 1:length(trans.t)-1
    tau_sim(i) = trans.tau(i) - trans.tau_prev;
    x_sim(i+1) = x_sim(i) + trans.h*(1/m*(tau_sim(i) + trans.k*(trans.ss(1) - x_sim(i))));
    error = error + (trans.x(i) - x_sim(i))^2; 
end

% MSE for the purpose of downscaling the problem to avoid extreme numbers
error = 1/length(trans.t) * error;
x     = x_sim;
end
