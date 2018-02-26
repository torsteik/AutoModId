function [m,k] = BFGS(trans)


% Declare and initialize
% - Iterator struct (BFGS):
%    - This struct contains all the information belonging to the
%      (outer-)loop, the main structure of the BFGS.
k = struct;
% - Design variables
k.low_bound = 0.1;            % Lower bound for inertia, m
k.m_delta   = 1e-8;           % Difference used in finite differnece
k.conv_tol  = 1e-5;           % Convergance tolerance (Stopping condition)
k.max_iter  = 50;

% - Core variables
k.m = NaN(1,k.max_iter);        % Inertia value (Optimization variable)

k.err     = NaN(1,k.max_iter);  % LS error (What is being minimized)
k.err_dot = NaN(1,k.max_iter);

k.p     = NaN(1,k.max_iter);    % Search direction
k.alpha = zeros(1,k.max_iter);  % Step length

% - Initialize search
k.m(1)     = k.low_bound;
k.err(1)   = eval_err(k.m(1),trans);

i = 1;
k.err_ddot = -1;
while k.err_ddot < 0
    k.err_delta  = eval_err(k.m(1) +  i*k.m_delta,trans);
    k.err_ddelta = eval_err(k.m(1) +2*i*k.m_delta,trans);

    k.err_dot(1) = (k.err_delta - k.err(1))/(i*k.m_delta);
    k.err_ddot   = (k.err_ddelta - 2*k.err_delta + k.err(1))/((i*k.m_delta)^2);
    
    i = i+1;
end
k.hessian_inv = 1/max(k.err_ddot,1e-2);

if trans.plot_bool == 1
    %% Plot        
    % Set up figure
    trans.e_fig = figure(302);
    clf(trans.e_fig,'reset')
    set(trans.e_fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
    
    % Plot of error
    trans.e_ax = subplot(2,1,1);
    xlabel(trans.e_ax,'Inertia, m');
    hold(trans.e_ax,'on');
    trans.e_ax.XGrid = 'on';
    trans.e_ax.YGrid = 'on';
    trans.e_ax.LineWidth = 1;
    
    trans.k_e_plot = gobjects(k.max_iter,1);
    trans.k_e_plot(1) = plot(trans.e_ax, k.m(1),k.err(1));
    trans.k_e_plot(1).Marker = 'o';
    trans.k_e_plot(1).MarkerSize = 6;
    trans.k_e_plot(1).Color = [1 0 0];
    
    % Plot of error derivative
    trans.ed_ax = subplot(2,1,2);
    xlabel(trans.ed_ax,'Inertia, m');
    hold(trans.ed_ax,'on');
    trans.ed_ax.XGrid = 'on';
    trans.ed_ax.YGrid = 'on';
    trans.ed_ax.LineWidth = 1;

%     trans.k_ed_plot = gobjects(k.max_iter,1);
%     trans.k_ed_plot(1) = plot(trans.ed_ax, k.m(1),k.err_dot(1));
%     trans.k_ed_plot(1).Marker = 'o';
%     trans.k_ed_plot(1).MarkerSize = 6;
%     trans.k_ed_plot(1).Color = [1 0 0];
    input(strcat('Loop coordinates: ts = ',num2str(trans.id),', k=',num2str(1),'\nPress ENTER to continue'))
end

% Start BFGS
k.k = 1;
while (k.k == 1 || abs(k.alpha(k.k-1)*k.p(k.k-1)) > k.conv_tol) &&  k.k < k.max_iter   % Additional condition is added in original. abs(k.err_dot(k.k))

    k.p(k.k)     = -k.hessian_inv*k.err_dot(k.k);
    k.alpha(k.k) = line_search(k,trans);            % Line search algorithm
    k.m(k.k + 1) = k.m(k.k) + k.alpha(k.k)*k.p(k.k);

    k.err(k.k+1)     = eval_err(k.m(k.k+1),trans);
    k.err_delta      = eval_err(k.m(k.k+1) + k.m_delta,trans);
    k.err_dot(k.k+1) = (k.err_delta - k.err(k.k+1))/k.m_delta;

    
    
    s_k   = k.m(k.k+1) - k.m(k.k);
    y_k   = k.err_dot(k.k+1) - k.err_dot(k.k);
    psi_k = 1/(y_k*s_k);
    k.hessian_inv  = (1-psi_k*s_k*y_k)*k.hessian_inv*(1-psi_k*y_k*s_k) + psi_k*s_k*s_k;
    k.hessian_inv2 = (k.m(k.k+1) - k.m(k.k))/(k.err_dot(k.k+1) - k.err_dot(k.k));

    k.k = k.k + 1;
    
    if trans.plot_bool == 1
        %% Plot
        % Recolor previous iterations
        trans.k_e_plot(k.k-1).Color     = [0 0 1];
        trans.k_ed_dot_plot(k.k-1).Color = [0 0 1];
        
        % Plot of error
        trans.k_e_plot(k.k) = plot(trans.e_ax, k.m(k.k),k.err(k.k));
        trans.k_e_plot(k.k).Marker = 'o';
        trans.k_e_plot(k.k).MarkerSize = 6;
        trans.k_e_plot(k.k).Color = [0 1 0];

        % Plot of error derivative
%         trans.k_ed_plot(k.k) = plot(trans.ed_ax, k.m(k.k),k.err_dot(k.k));
%         trans.k_ed_plot(k.k).Marker = 'o';
%         trans.k_ed_plot(k.k).MarkerSize = 6;
%         trans.k_ed_plot(k.k).Color = [0 1 0];
        input(strcat('Loop coordinates: ts = ',num2str(trans.id),', k=',num2str(k.k),'\nPress ENTER to continue'))
    end
end

m = k.m(k.k);
end