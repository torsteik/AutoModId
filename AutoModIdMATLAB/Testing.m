% SCRIPT #2
load('Experiment.mat');
%% Variable and parameter declarations
N = length(exp.time);
h = exp.h;

n = exp.type + 1;

% Gaussian noise
x_noisy = zeros(1,N);

noise  = zeros(1,N);
random = rand(2,N);
stdev  = (n == 1)*1e-6 + (n == 2)*1e-7;

% Key variables
x_f     = zeros(1,N);
nu_f    = zeros(1,N);
delta_f = zeros(1,N);

R = zeros(1,N);

% Filter factors and variance thresholds
% NOE MÅ GJØRES MED n HER
lambda_1 = 0.15*h;
lambda_2 = 0.5*h;
lambda_3 = 0.1*h;

R_ss = 5000;
R_ts = 8000;

% Steady state information variables
ss_flag       = zeros(1,N);                  % Logs steady state: TRANSIENT = 0, STEADY_STATE = 1
ss_timestamps = ones(2,length(exp.tau_seq)); % Logs time of steady state: [ENTRANCE, EXIT]

%ss_flag(1)  = 1;
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
            continue;
        
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
        if max(exp.tau(:,i) ~= exp.tau(:,i-1))
            ss_flag(i) = 0;
            ss_timestamps(2,tau_seq_itr) = i;
            tau_seq_itr = tau_seq_itr +1;
            
            if exp.tau(n,i) == exp.tau(n,i-1)
                ts_bool = 1;
            end
        else
            ss_flag(i) = 1;
        end
    end
    
    if tau_seq_itr == length(exp.tau_seq) +1
        break;
    end
end


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
ylabel(R_ax,'$R$')
R_R_plot.LineWidth = 1;
%R_ax.YLim = [0 1];

yyaxis(R_ax,'left')
R_x_plot   = plot(R_ax, exp.time, exp.x(n,:).*(r2d^(n-1)));
R_x_f_plot = plot(R_ax, exp.time, x_f.*(r2d^(n-1)),'--k');
ylabel(R_ax,strcat((n==1)*'SOG [m/s]',(n==2)*'ROT [deg/s]'))
R_x_plot.LineWidth   = 1;
R_x_f_plot.LineWidth = 1;
R_ax.YLim = [0 max(exp.x(n,:))*(r2d^(n-1))];


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
%R_ax.YLim = [0 1];

yyaxis(SS_ax,'left')
SS_area_plot = area(SS_ax, exp.time, ss_flag.*(max(exp.x(n,:).*(r2d^(n-1)))+1));
SS_line_plot = plot(SS_ax, exp.time, ss_flag.*(max(exp.x(n,:).*(r2d^(n-1)))+1));
SS_x_plot  = plot(SS_ax, exp.time, exp.x(n,:).*(r2d^(n-1)));
ylabel(SS_ax,strcat((n==1)*'SOG [m/s]',(n==2)*'ROT [deg/s]'))
SS_ax.YLim = [0 max(exp.x(n,:))*(r2d^(n-1))];
SS_x_plot.LineWidth  = 1;
SS_x_plot.LineStyle  = '-';
SS_area_plot.LineWidth = 1;
SS_area_plot.LineStyle = '-';
SS_area_plot.FaceColor = [0.75 0.75 0.75];
SS_area_plot.FaceAlpha = 0.3;
SS_line_plot.LineStyle = '-';
SS_line_plot.Color = 'k';


%% NOTES
% Maybe include option to only perform the calculation of R. For example
% use 'input()' to ask and then put an 'if()' after the R calculation.