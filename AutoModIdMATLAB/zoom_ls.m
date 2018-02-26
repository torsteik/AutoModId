function alpha = zoom_ls(lo, hi, k, i, trans)
% Interpolate to find step lenght between low.alpha and high.alpha. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lo.alpha   = lo.alpha(lo.i);
lo.err     = lo.err(lo.i);
lo.err_dot = lo.err_dot(lo.i);

hi.alpha   = hi.alpha(hi.i);
hi.err     = hi.err(hi.i);
hi.err_dot = hi.err_dot(hi.i);

% Iterator struct (Zoom):
% - This struct contains all the information belonging to the
%   (inner-)loop, the main structure of the zoom.
j = struct;

j.max_iter = 50;

j.alpha   = NaN(1,j.max_iter);
j.err     = NaN(1,j.max_iter);
j.err_dot = NaN(1,j.max_iter);

if trans.plot_bool == 1
    %% Plot
    % Plot of error
    trans.j_e_plot = gobjects(j.max_iter+1,1);
    % - Plot lo
    trans.j_e_plot(1) = plot(trans.e_ax, k.m + lo.alpha*k.p, lo.err);
    trans.j_e_plot(1).Marker = 'x';
    trans.j_e_plot(1).MarkerSize = 8;
    trans.j_e_plot(1).Color = [1 0.75 0];
    % - Plot hi
    trans.j_e_plot(2) = plot(trans.e_ax, k.m + hi.alpha*k.p, hi.err);
    trans.j_e_plot(2).Marker = 'x';
    trans.j_e_plot(2).MarkerSize = 8;
    trans.j_e_plot(2).Color = [1 0 1];
    
    % Plot of error derivative
    trans.j_ed_plot = gobjects(j.max_iter+1,1);
    % - Plot lo
    trans.j_ed_plot(1) = plot(trans.ed_ax, k.m + lo.alpha*k.p, lo.err_dot);
    trans.j_ed_plot(1).Marker = 'x';
    trans.j_ed_plot(1).MarkerSize = 8;
    trans.j_ed_plot(1).Color = [1 0.75 0];
    % - Plot hi
    trans.j_ed_plot(2) = plot(trans.ed_ax, k.m + hi.alpha*k.p, hi.err_dot);
    trans.j_ed_plot(2).Marker = 'x';
    trans.j_ed_plot(2).MarkerSize = 8;
    trans.j_ed_plot(2).Color = [1 0 1];
    input(strcat('Loop coordinates: k=',num2str(k.k),', i=',num2str(i.i),', and j=',num2str(1),...
                 '\nPress ENTER to continue'))
end

% Start zoom
j.j = 1;
while (j.j == 1 || abs(j.alpha(j.j-1)*k.p) > k.conv_tol) && j.j <= j.max_iter
    % Quadratic interpolation - Assuming it's okay to set a_i = a_hi and a_i-1 = a_lo
%     if lo.alpha < hi.alpha
%         j.alpha(j.j) = -lo.err_dot*k.p*hi.alpha^2/(2*(hi.err-lo.err-lo.err_dot*k.p*hi.alpha));
%     else
%         j.alpha(j.j) = -hi.err_dot*k.p*lo.alpha^2/(2*(lo.err-hi.err-hi.err_dot*k.p*lo.alpha));
%     end
    j.alpha(j.j) = 0.5*(lo.alpha + hi.alpha);
    % Evaluate objective
    % and estimate derivative of objective in search direction
    j.err(j.j) = eval_err(k.m + j.alpha(j.j)*k.p, trans);
    j.err_delta    = eval_err((k.m + k.m_delta) + j.alpha(j.j)*k.p, trans);
    j.err_dot(j.j) = (j.err_delta - j.err(j.j))/k.m_delta;
    
    if trans.plot_bool == 1
        %% Plot
        
        input('Press ENTER to continue')
    end
    
    if trans.plot_bool == 1
        %% Plot
        % Plot of error
        trans.j_e_plot(j.j+1) = plot(trans.e_ax, k.m + j.alpha(j.j)*k.p, j.err(j.j));
        trans.j_e_plot(j.j+1).Marker = 'x';
        trans.j_e_plot(j.j+1).MarkerSize = 8;
        trans.j_e_plot(j.j+1).Color = [0 0 0];

        % Plot of error derivative
        trans.j_ed_plot(j.j+1) = plot(trans.ed_ax, k.m + j.alpha(j.j)*k.p, j.err_dot(j.j));
        trans.j_ed_plot(j.j+1).Marker = 'x';
        trans.j_ed_plot(j.j+1).MarkerSize = 8;
        trans.j_ed_plot(j.j+1).Color = [0 0 0];
        
        % Plog transient
        r2d = 180/pi;
        [~,x_sim] = eval_err(k.m + j.alpha(j.j)*k.p,trans);
        trans.ts_plot = plot(trans.ts_ax, trans.t.*trans.h, x_sim*(r2d^(trans.n-1)));
        
        input(strcat('Loop coordinates: k=',num2str(k.k),', i=',num2str(i.i),', and j=',num2str(j.j),...
                     '\nPress ENTER to continue'))
    end

    % Check sufficient decrease condition
    if j.err(j.j) > k.err + i.c_1*j.alpha(j.j)*k.err_dot*k.p || j.err(j.j) >= lo.err
        hi.alpha = j.alpha(j.j);
        
        if trans.plot_bool == 1
            %% Plot
            % Change color of point to hi-color
            trans.j_e_plot(j.j+1).Color = [1 0 1];
            input('Press ENTER to continue')
        end
    else
        % Check curvature condition
        if abs(j.err_dot(j.j)*k.p) <= -i.c_2*k.err_dot*k.p
            alpha = j.alpha(j.j);
            break;
        end
        
        % Similar interpretation to: 
        % Stepped to far but within sufficient decrease region
        if j.err_dot(j.j)*k.p*(hi.alpha - lo.alpha) >= 0
            hi.alpha = lo.alpha;
        end
        
        lo.alpha = j.alpha(j.j);
        if trans.plot_bool == 1
            %% Plot
            % Change color of point to lo-color
            trans.j_e_plot(j.j+1).Color = [1 0.75 0];
            input('Press ENTER to continue')
        end
    end
    
    j.j = j.j + 1;
end


if j.j > j.max_iter || exist('alpha') ~= 1
    alpha = j.alpha(j.j-1);
end
end