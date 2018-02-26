% SCRIPT #2
n = input('Do you want to load:\n1 - Motor experiment\n2 - Rudder experiment\nAnswer: ');
if n == 1
    load('ExperimentMotor.mat');
elseif n == 2
    load('ExperimentRudder.mat');
else
    'Invalid input'
end
%% Variable and parameter declarations
N = length(exp.time);
h = exp.h;

n = exp.type + 1;

% Gaussian noise
x_noisy = zeros(1,N);

noise  = zeros(1,N);
random = rand(2,N);
stdev  = 0;%(n == 1)*1e-6 + (n == 2)*1e-7;

% Key variables
x_f     = zeros(1,N);
nu_f    = zeros(1,N);
delta_f = zeros(1,N);

R = zeros(1,N);

% Filter factors and variance thresholds
if n == 1
    lambda_1 = 0.2*h; % 0.2 og 0.09 Ser også bra ut
    lambda_2 = 0.5*h;
    lambda_3 = 0.075*h;   
    
    R_ss = 700;%7000(lavere med mer tid) og 2000;
    R_ts = 2000;%10000 og 4000;
else
    lambda_1 = 0.3*h;
    lambda_2 = 0.75*h;
    lambda_3 = 0.5*h;
    
    R_ss = 2000;
    R_ts = 6000;
end


% Steady state information variables
ss_flag       = zeros(1,N);                  % Logs steady state: TRANSIENT = 0, STEADY_STATE = 1
ss_timestamps = ones(2,length(exp.tau_seq)); % Logs time of steady state: [ENTRANCE, EXIT]

ts_bool  = 1; % Used to ensure transient gets going before SSID starts
tau_seq_itr = 1; % Keeping track of how far into the experiment we are
%% Simulation of experiment
for i = 2:N
    % Add noise to state measurements from experiment
    noise(i) = stdev.*(sqrt(-2*log(random(1,i)))*sin(2*pi*random(2,i)))';
    x_noisy(i) = exp.x(n,i) + noise(i);
    
    % Filter state
    x_f(i) = lambda_1*x_noisy(i) + (1 - lambda_1)*x_f(i-1);
    
    % Variance estimates
    nu_f(i)    = lambda_2*(x_noisy(i)-x_f(i-1))^2     + (1-lambda_2)*nu_f(i-1);
    delta_f(i) = lambda_3*(x_noisy(i)-x_noisy(i-1))^2 + (1-lambda_3)*delta_f(i-1);
    
    R(i) = ((2-lambda_1)*nu_f(i)) / delta_f(i);
    
    % Switch
    % - Case 1: Transient
    if ss_flag(i-1) == 0
        if R(i) < R_ts && ts_bool == 0
            %Do nothing
        
        else
            ts_bool = 1;
            if R(i) < R_ss
                ts_bool = 0;
                ss_flag(i) = 1;
                ss_timestamps(1,tau_seq_itr) = i;
            end
        end
    end
    
    % - Case 2: Steady state
    if ss_flag(i-1) == 1
        if any(exp.tau(:,i) ~= exp.tau(:,i-1))
            ss_flag(i) = 0;
            % Commented out in order to handle poor data
            % ss_timestamps(2,tau_seq_itr) = i;
            % tau_seq_itr = tau_seq_itr +1;
            
            if exp.tau(n,i) == exp.tau(n,i-1)
                ts_bool = 1;  % Should do something smarter here maybe
            end
        else
            ss_flag(i) = 1;
        end
    end
 
    % Replacement for commented out lines:
    % Skips to next input pair in tau_seq when new input is applied,
    % regardless of whether or not steady state was achieved or not.
    if any(exp.tau(:,i) ~= exp.tau(:,i-1))
        ss_timestamps(2,tau_seq_itr) = i;
        tau_seq_itr = tau_seq_itr +1;
        if ~(exp.tau(n,i) == exp.tau(n,i-1))
                ts_bool = 0;
        end
    end
    
    if tau_seq_itr == length(exp.tau_seq) +1
        break;
    end
end

ss_timestamps(2,end) = length(exp.time);

%% Plot
set(groot, 'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

r2d = 180/pi;
% Plot #1: R-statistic
R_fig = figure(201);
clf(R_fig,'reset')
set(R_fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);

R_ax = axes;
xlabel(R_ax,'time [s]');
hold(R_ax,'on');
R_ax.XGrid = 'on';
R_ax.YGrid = 'on';
R_ax.LineWidth = 1;
R_ax.XLim = [exp.time(1) exp.time(end)];

yyaxis(R_ax,'right')
R_R_plot = plot(R_ax, exp.time, R);
R_nu_plot = plot(R_ax, exp.time, nu_f);
R_delta_plot = plot(R_ax, exp.time, delta_f);
ylabel(R_ax,'$R$')
R_R_plot.LineWidth = 1;
R_nu_plot.LineWidth = 1;
R_nu_plot.LineStyle = '-';
R_nu_plot.Color = [0 0.5977 0.1992];
R_detla_plot.LineWidth = 1;
R_delta_plot.LineStyle = '-';
R_delta_plot.Color = [1 0.5977 0];
%R_ax.YLim = [0 0.5e5];

yyaxis(R_ax,'left')
R_x_plot   = plot(R_ax, exp.time, exp.x(n,:).*(r2d^(n-1)));
R_x_f_plot = plot(R_ax, exp.time, x_f.*(r2d^(n-1)),'--k');
R_tau_plot = plot(R_ax, exp.time, exp.tau(n,:),'-g');
ylabel(R_ax,strcat((n==1)*'SOG [m/s]',(n==2)*'ROT [deg/s]'))
R_x_plot.LineWidth   = 1;
R_x_f_plot.LineWidth = 1;
R_tau_plot.LineWidth = 1;
%R_ax.YLim = [0 2];
R_ax.YLim = [0 max(exp.x(n,:))*(r2d^(n-1))];

legend('$x$','$x_f$','$R$','$\nu$','$\delta$')


% Plot #2: SSID result
SS_fig = figure(202);
clf(SS_fig,'reset')
set(SS_fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);

SS_ax = axes;
xlabel(SS_ax,'time [s]');
hold(SS_ax,'on');
SS_ax.XGrid = 'off';
SS_ax.YGrid = 'off';
SS_ax.LineWidth = 1;
SS_ax.XLim = [exp.time(1) exp.time(end)];

yyaxis(SS_ax,'right')
SS_R_plot = plot(SS_ax, exp.time, R);
ylabel(SS_ax,'$R$')
SS_R_plot.LineWidth = 1;
%SS_ax.YLim = [0 2.5e4];

yyaxis(SS_ax,'left')
SS_area_plot = area(SS_ax, exp.time, ss_flag.*(max(exp.x(n,:).*(r2d^(n-1)))+1));
SS_line_plot = plot(SS_ax, exp.time, ss_flag.*(max(exp.x(n,:).*(r2d^(n-1)))+1));
SS_x_plot    = plot(SS_ax, exp.time, exp.x(n,:).*(r2d^(n-1)));
SS_tau_plot  = plot(SS_ax, exp.time, exp.tau(n,:),'-g');
ylabel(SS_ax,strcat((n==1)*'SOG [m/s]',(n==2)*'ROT [deg/s]'))
%SS_ax.YLim = [0 2];
SS_ax.YLim = [0 max(exp.x(n,:))*(r2d^(n-1))];
SS_x_plot.LineWidth  = 1;
SS_x_plot.LineStyle  = '-';
SS_area_plot.LineWidth = 1;
SS_area_plot.LineStyle = '-';
SS_area_plot.FaceColor = [0.75 0.75 0.75];
SS_area_plot.FaceAlpha = 0.3;
SS_line_plot.LineStyle = '-';
SS_line_plot.Color = 'k';

%% Save and exit
if n == 1
    save('ActualExperimentDataMotor.mat','exp','ss_timestamps');
else
    save('ActualExperimentDataRudder.mat','exp','ss_timestamps');
end
