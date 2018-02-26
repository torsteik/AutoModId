function alpha = line_search(k,trans)
% Find step length satisfying wolfe conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k.m         = k.m(k.k);
k.p         = k.p(k.k);
k.alpha     = k.alpha(k.k);  % Just a zero at this point
k.err       = k.err(k.k);
k.err_dot = k.err_dot(k.k);

% 'k.err_dot' represents the derivative of the obj. func. with respect to
% the opt. var., 'm'. But from here on derivatives are considered with
% respect to the step length, 'alpha'.
% k.err_delta = eval_err(k.m + (k.alpha + k.m_delta)*k.p,trans);
% k.err_dot   = (k.err_delta - k.err)/k.m_delta;

% Iterator struct (Line search):
% - This struct contains all the information belonging to the
%   (middle-)loop, the main structure of the line search.
i = struct;

i.max_iter  = 50;

i.alpha   = NaN(1,i.max_iter);
i.err     = NaN(1,i.max_iter);
i.err_dot = NaN(1,i.max_iter);

i.alpha(1)   = k.alpha;
i.err(1)     = k.err;
i.err_dot(1) = k.err_dot;

% Design variables
i.alpha_min = 0;
i.alpha_max = 10;

i.c_1 = 0.1;      % Wolfe condition parameter
i.c_2 = 0.9;      % ----------- " -----------

i.alpha(2) = 1;

if trans.plot_bool == 1
    %% Plot
    % Plot of error
    trans.i_e_plot = gobjects(i.max_iter,1);
    trans.i_e_plot(1) = plot(trans.e_ax, k.m,i.err(1));
    trans.i_e_plot(1).Marker = '.';
    trans.i_e_plot(1).MarkerSize = 8;
    trans.i_e_plot(1).Color = [1 0 0];

    % Plot of error derivative
    trans.i_ed_plot = gobjects(i.max_iter,1);
    trans.i_ed_plot(1) = plot(trans.ed_ax, k.m,i.err_dot(1));
    trans.i_ed_plot(1).Marker = '.';
    trans.i_ed_plot(1).MarkerSize = 8;
    trans.i_ed_plot(1).Color = [1 0 0];

    input(strcat('Loop coordinates: k=',num2str(k.k),', i=',num2str(1),'\nPress ENTER to continue'))
end

% Start line search
i.i = 2;
while (i.i == 2 || abs(i.alpha(i.i-1)*k.p) > k.conv_tol) && i.i <= i.max_iter  % These conditions should not be necessary.
    % Evaluate objective at initial step
    % and estimate derivative of objective in search direction
    i.err(i.i)     = eval_err(k.m + i.alpha(i.i)*k.p,trans);
    i.err_delta    = eval_err((k.m + k.m_delta) + i.alpha(i.i)*k.p,trans);
    i.err_dot(i.i) = (i.err_delta - i.err(i.i))/k.m_delta;
    
    if trans.plot_bool == 1
        %% Plot
        % Recolor previous iterations
        trans.i_e_plot(i.i-1).Color      = [0 0 1];
        trans.i_ed_dot_plot(i.i-1).Color = [0 0 1];

        % Plot of error
        trans.i_e_plot(i.i) = plot(trans.e_ax, k.m + i.alpha(i.i)*k.p, i.err(i.i));
        trans.i_e_plot(i.i).Marker = '.';
        trans.i_e_plot(i.i).MarkerSize = 8;
        trans.i_e_plot(i.i).Color = [0 0.75 0];

        % Plot of error derivative
        trans.i_ed_plot(i.i) = plot(trans.ed_ax, k.m + i.alpha(i.i)*k.p, i.err_dot(i.i));
        trans.i_ed_plot(i.i).Marker = '.';
        trans.i_ed_plot(i.i).MarkerSize = 8;
        trans.i_ed_plot(i.i).Color = [0 0.75 0];
        
        % Plot of transient
        r2d = 180/pi;
        [~,x_sim] = eval_err(k.m + i.alpha(i.i)*k.p,trans);
        trans.ts_plot = plot(trans.ts_ax, trans.t.*trans.h, x_sim*(r2d^(trans.n-1)));
        
        input(strcat('Loop coordinates: k=',num2str(k.k),', i=',num2str(i.i),'\nPress ENTER to continue'))
    end
    
    % Check sufficient decrease condition
    if i.err(i.i) > k.err + i.c_1*i.alpha(i.i)*k.err_dot*k.p || (i.i > 2 && i.err(i.i) >= i.err(i.i-1))
        % Declare 'lo' to be the previous iteration of 'i' and 'hi' the current
        lo   = i; 
        lo.i = i.i - 1;
        
        hi = i;
        % Zoom
        alpha = zoom_ls(lo, hi, k, i, trans);
        break;
    end

    % Check curvature condition
    if abs(i.err_dot(i.i)*k.p) <= -i.c_2*k.err_dot*k.p
        % A solution has been identified
        alpha = i.alpha(i.i);
        break;
    end

    % Stepped to far but within sufficient decrease region
    if i.err_dot(i.i)*k.p >= 0
        % Declare 'hi' to be the previous iteration of 'i' and 'lo' the current
        lo = i;
        
        hi   = i; 
        hi.i = i.i - 1;
        % Zoom
        alpha = zoom_ls(lo, hi, k, i, trans);
        break;
    end

    % Last possibility, insufficient step length
    % Extrapolate
    alpha = i.alpha(i.i);
    i.alpha(i.i+1) = 5*i.alpha(i.i);
    
    i.i = i.i + 1;
end

if i.i > i.max_iter
    alpha = i.alpha(i.i-1);
end

end
