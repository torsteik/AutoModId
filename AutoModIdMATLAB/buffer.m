% SCRIPT #3
n = input('\nDo you want to load:\n1 - Motor experiment\n2 - Rudder experiment\nAnswer: ');
if n == 1
    load('ActualExperimentDataMotor.mat');
elseif n == 2
    load('ActualExperimentDataRudder.mat');
else
    'Invalid input'
end
n = exp.type +1;
%% Damping measurements
% Damping is found directly as the input
sigma_r = exp.tau_seq(2,:);
sigma_r(ismember(ss_timestamps(1,:),1)) = []; % Remove the invalid elements

% Get steady state averages
valid_timestamps = ss_timestamps(:,~ismember(ss_timestamps(1,:),1));
x_ss_avg = zeros(2,length(sigma_r));

for i = 1:length(sigma_r)
    x_ss_avg(:,i) = mean((exp.x(:,valid_timestamps(1,i):valid_timestamps(2,i))'));
end

%% Inertia measurements

plot_bool   = input('\nTurn on intermediate plots?\n1 - Yes\n2 - No\nAnswer: ');
first_trans = input('\nSkip to specific transient?\nTransient number(1 = No): ');


% Go through all identified transients
for trans_itr = first_trans:length(valid_timestamps)-1
    % Skip invalid data and steps in secondary input
    if any(ismember(exp.secondary_steps,trans_itr)) || (ss_timestamps(1,trans_itr+1) == 1)
        continue;
    end
    
    %% Set up
    % Create struct to store info about the transient
    trans     = struct;
    trans.h   = exp.h;
    trans.t   = ss_timestamps(2,trans_itr):ss_timestamps(1,trans_itr+1);
    trans.x   = exp.x(n,trans.t);
    trans.tau = exp.tau(n,trans.t);
    trans.ss  = [x_ss_avg(n,trans_itr);
                x_ss_avg(n,trans_itr+1)];
    trans.tau_prev = exp.tau_seq(n,trans_itr);
    trans.k   = (sigma_r(trans_itr+1) - sigma_r(trans_itr))/(x_ss_avg(n,trans_itr+1) - x_ss_avg(n,trans_itr)); % Linearized damping term
    
    if plot_bool == 1
        %% Plot
        r2d = 180/pi;
        
        ts_fig = figure(301);
        clf(ts_fig,'reset')
        set(ts_fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
        
        ts_ax = axes;
        xlabel(ts_ax,'time [s]');
        hold(ts_ax,'on');
        ts_ax.XGrid = 'on';
        ts_ax.YGrid = 'on';
        ts_ax.LineWidth = 1;
        
        ts_plot = plot(ts_ax, trans.t.*trans.h, trans.x.*(r2d^(n-1)));
        ts_ax.XLim = [trans.t(1) trans.t(end)].*trans.h;
        
        input('Press ENTER to continue')
    end
    %% BFGS
    % Declare and initialize
    % - Iterator struct (BFGS):
    %    - This struct contains all the information belonging to the
    %      (outer-)loop, the main structure of the BFGS.
    k = struct;
    % - Design variables
    k.low_bound = 0.1;            % Lower bound for inertia, m
    k.m_delta   = 1e-8;           % Difference used in finite differnece
    k.conv_tol  = 1e-8;           % Convergance tolerance (Stopping condition)
    k.max_iter  = 50;
    
    % - Core variables
    k.m = NaN(1,k.max_iter);        % Inertia value (Optimization variable)
    
    k.err = NaN(1,k.max_iter);      % LS error (What is being minimized)
    k.err_dot = NaN(1,k.max_iter);
    
    k.p     = NaN(1,k.max_iter);    % Search direction
    k.alpha = NaN(1,k.max_iter);    % Step length
    
    % - Initialize search
    k.m(1)     = k.low_bound;
    k.alpha(1) = 0;
    k.err(1)   = eval_err(m(1),trans);
    
    k.err_delta  = eval_err(m(1) +  k.m_delta,trans);
    k.err_ddelta = eval_err(m(1) +2*k.m_delta,trans);
    
    k.err_dot(1) = (k.err_delta - k.err(1))/k.m_delta;
    k.err_ddot   = (k.err_ddelta - 2*k.err_delta + k.err(1))/(k.m_delta^2);
    
    k.hessian_inv = 1/k.err_ddot;
    
    if plot_bool == 1
        %% Plot        
        % Set up figure
        err_fig = figure(302);
        clf(err_fig,'reset')
        set(err_fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
        
        err_ax = axes;
        xlabel(err_ax,'Inertia, m');
        hold(err_ax,'on');
        err_ax.XGrid = 'on';
        err_ax.YGrid = 'on';
        err_ax.LineWidth = 1;
        
        init_plot = plot(err_ax, k.m(1),k.err(1));
        init_plot.Marker = 'o';
        init_plot.MarkerSize = 6;
        init_plot.Color = [1 0 0];
    end
    
    % Start BFGS
    k.k = 1;
    while abs(k.err_dot(k.k)) > k.conv_tol &&  k.k < k.max_iter   % Additional condition is added in original
        
        % Compute search direction
        k.p(k.k) = -k.hessian_inv*k.err_dot(k.k);
        
        %% Line search
        % Find step length using a line search algorithm
        
        % Decleration
        % - Iterator struct (Line search):
        %    - This struct contains all the information belonging to the
        %      (middle-)loop, the main structure of the line search.
        i = struct;
        
        i.max_iter  = 50;
        i.alpha   = NaN(1,i.max_iter);
        i.err     = NaN(1,i.max_iter);
        i.err_dot = NaN(1,i.max_iter);
        
        % Design variables
        i.alpha_min = 0;
        i.alpha_max = 10;
        
        i.c_1 = 0.1;      % Wolfe condition parameter
        i.c_2 = 0.9;      % ----------- " -----------
         
        i.alpha(1) = 1;
        
        % Start line search
        i.i = 1;
        while abs(i.err_dot(i.i)) > k.conv_tol && i.i < i.max_iter  % These conditions should not be necessary.
            % Evaluate objective at initial step
            i.err(i.i) = eval_err(m(k.k) + i.alpha(i.i)*k.p(k.k),trans);
            
            % Check sufficient decrease condition
            if i.err(i.i) > k.err(k.k) + i.c_1*i.alpha(i.i)*k.err_dot(k.k) || (i.err(i.i) >= i.err-1)
                % DO SHIT
            end
            
            % Estimate derivative of objective in search direction
            i.err_delta  = eval_err(m(k.k) + (i.alpha(i.i)+k.m_delta)*k.p(k.k),trans);
            i.err_dot(i.i) = (i.err_delta - i.err(i.i))/k.m_delta;
            
            % Check curvature condition
            if abs(i.err_dot(i.i)) <= -i.c_2*k.err_dot(k.k)
                % DO OTHER SHIT
            end
            
            % Stepped to far but within sufficient decrease region
            if i.err_dot(i.i) >= 0
                % DO SIMILAR SHIT TO THE FIRST ONE
            end
            
            % Last possibility, insufficient step length
            % EXTRAPOLATE
            
            i.i = i.i + 1;
        end
        
         k.k = k.k + 1;
    end
end




% 
% 
% figure(1)
% hold on
% grid on
% scatter3(x_ss(2,cv_ind == 1), x_ss(1,cv_ind == 1), sigma(cv_ind == 1));
% scatter3(x_ss(2,cv_ind == 2), x_ss(1,cv_ind == 2), sigma(cv_ind == 2));
% scatter3(x_ss(2,cv_ind == 3), x_ss(1,cv_ind == 3), sigma(cv_ind == 3));
% scatter3(x_ss(2,cv_ind == 4), x_ss(1,cv_ind == 4), sigma(cv_ind == 4));
% scatter3(x_ss(2,cv_ind == 5), x_ss(1,cv_ind == 5), sigma(cv_ind == 5));
% scatter3(x_ss(2,cv_ind == 6), x_ss(1,cv_ind == 6), sigma(cv_ind == 6));
% scatter3(x_ss(2,cv_ind == 7), x_ss(1,cv_ind == 7), sigma(cv_ind == 7));
% scatter3(x_ss(2,cv_ind == 8), x_ss(1,cv_ind == 8), sigma(cv_ind == 8));
% scatter3(x_ss(2,cv_ind == 9), x_ss(1,cv_ind == 9), sigma(cv_ind == 9));
% scatter3(x_ss(2,cv_ind == 10), x_ss(1,cv_ind == 10), sigma(cv_ind == 10));




