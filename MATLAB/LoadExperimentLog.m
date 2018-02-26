% SCRIPT #1
%% Read logs from experiment.
load = input('Load experiment?\n 1 - Yes\n 2 - No\nAnswer: ');
if load == 1
    clear;
    addpath(genpath('ExperimentData'));
    
    % Experiments:
    %  - 'MR_log_surge.log' + 'SeaPath_log_surge.log' 
    %    => motor step experiment
    %  - 'MR_log_rudder.log' + ('SeaPath_log_surge_1.log' and 'SeaPath_log_surge_2.log')
    %    => rudder step experiment, where SeaPath_log_surge_2.log is the best
    experiment = input('Choose experiment\n 1 - surge\n 2 - rudder1\n 3 - rudder2\nAnswer: ');

    MR = struct;
    SP = struct;
    if experiment == 1
        MR.log = load('MR_log_surge.log');
        SP.log = load('SeaPath_log_surge.log');
    elseif experiment == 2
        MR.log = load('MR_log_rudder.log');
        SP.log = load('SeaPath_log_rudder_1.log');
    elseif experiment == 3
        MR.log = load('MR_log_rudder.log');
        SP.log = load('SeaPath_log_rudder_2.log');
    else
        'Invalid input'
        return;
    end

    % Log indices
    MR.time_ind = 1;
    SP.time_ind = [3 4];

    MR.tau_delta_ind = 4;
    MR.tau_m_ind     = 6;

    SP.r_ind = 17;

    SP.north_vel_ind = 9;
    SP.east_vel_ind  = 10;

    % Experiment parameters
    h = 0.01;
end

%% Interpret log
% Time - SeaPath
SP.time_abs = sum(SP.log(:,SP.time_ind),2);
SP.time_abs_start = SP.time_abs(1);
SP.time_abs_end   = SP.time_abs(end); 

SP.time = SP.time_abs - SP.time_abs_start;

% Time - MaritimeRobotics
% NOTE: Here approximately 0.04 seconds of inaccuracy is added, but this 
%       should be insignificant.
MR.time_start = find(MR.log(:,MR.time_ind) >= SP.time_abs_start,1,'first');
MR.time_end   = find(MR.log(:,MR.time_ind) < SP.time_abs_end,  1,'last');

MR.time_abs = MR.log(MR.time_start:MR.time_end, MR.time_ind);

MR.time = MR.time_abs - SP.time_abs_start;

% Input
MR.tau_m     = MR.log(MR.time_start:MR.time_end, MR.tau_m_ind);
MR.tau_delta = MR.log(MR.time_start:MR.time_end, MR.tau_delta_ind);

MR.tau_delta = MR.tau_delta./40;%(MR.tau_delta - 4)./36;  % Remove rudder offset and downscale

% State
SP.U = sqrt( sum( SP.log(:, [SP.north_vel_ind SP.east_vel_ind]).^2, 2) );
SP.r = SP.log(:, SP.r_ind);

%% Convert and store data in struct
exp = struct;

exp.h = h;
exp.type = (experiment > 1);    % 0 - motor steps, 1 - rudder steps

% Store time and state
exp.time = SP.time; 
exp.x    = [SP.U'; SP.r'];

% Upsample and store input
tau_m     = zeros(length(SP.time),1); 
tau_delta = zeros(length(SP.time),1);

j = 1;
for i = 2:length(MR.time)
    for j = j:length(SP.time)
        tau_m(j)     = MR.tau_m(i-1);
        tau_delta(j) = MR.tau_delta(i-1);
        if SP.time(j) >= MR.time(i)
            j = j + 1;
            break
        end
    end
end

exp.tau = [tau_m'; tau_delta'];
exp.tau_seq = [];                   % Input sequence. Initiated at the end.

%% Plot
set(groot, 'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

edit_exp = 1;
while edit_exp == 1
    % Plot experiment
    U_fig = figure(100+exp.type);
    clf(U_fig,'reset')
    set(U_fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);

    % SOG
    U_ax = subplot(2,1,1);
    xlabel(U_ax,'time [s]');
    U_ax.XGrid = 'on';
    U_ax.YGrid = 'on';
    U_ax.LineWidth = 1;
    U_ax.XLim = [exp.time(1) exp.time(end)];

    yyaxis(U_ax,'right')
    U_tau_plot = plot(U_ax, exp.time, exp.tau(1,:));
    ylabel(U_ax,'$\tau_m$')
    U_tau_plot.LineWidth = 1;
    U_ax.YLim = [0 1];

    yyaxis(U_ax,'left')
    U_x_plot = plot(U_ax,   exp.time, exp.x(1,:));
    ylabel(U_ax,'SOG [m/s]')
    U_x_plot.LineWidth = 1;
    U_ax.YLim = [0 max(exp.x(1,:))];


    % ROT
    r2d = 180/pi;

    r_ax = subplot(2,1,2);
    xlabel(r_ax,'time [s]');
    r_ax.XGrid = 'on';
    r_ax.YGrid = 'on';
    r_ax.LineWidth = 1;
    r_ax.XLim = [exp.time(1) exp.time(end)];

    yyaxis(r_ax,'right')
    r_tau_plot = plot(r_ax, exp.time, exp.tau(2,:));
    ylabel(r_ax,'$\tau_{\delta}$')
    r_tau_plot.LineWidth = 1;
    r_ax.YLim = [0 1];

    yyaxis(r_ax,'left')
    r_x_plot = plot(r_ax,   exp.time, exp.x(2,:)*r2d);
    ylabel(r_ax,'ROT [deg/s]')
    r_x_plot.LineWidth = 1;
    r_ax.YLim = [0 max(exp.x(2,:))*r2d];
    
    % Remove intervals from experiment
    edit_exp = input('Remove interval from experiment?\n 1 - Yes\n 2 - No\nAnswer: ');
    if edit_exp == 1
        start = input('Enter starting point of interval. [seconds]\nStart: ')/h;
        stop  = input('Enter stopping point of interval. [seconds]\nStop: ')/h;

        start = max(start,1);
        stop  = min(stop,length(exp.time));
        
        exp.time((end-(stop-start)):end) = [];
        exp.x(:,start:stop)   = [];
        exp.tau(:,start:stop) = [];
    end
end

%% Identify and store input sequence
% Edge detector - PogChamp
edge_m     = abs(exp.tau(1,2:end) - exp.tau(1,1:(end-1)));
edge_delta = abs(exp.tau(2,2:end) - exp.tau(2,1:(end-1))); 

edge = [(edge_m + edge_delta) 1];

exp.tau_seq = exp.tau(:, edge > 0);

% Find steps in secondary input (e.g. rudder step for motor experiment)
n = exp.type;

buffer = [1 find(edge_m*(n==1)+edge_delta*(n==0))];
exp.secondary_steps = zeros(1,length(buffer)-1);

tau_seq_itr = 0;
for i = 2:length(buffer)
    tau_seq_itr = 1 + tau_seq_itr + length(find((n==0)*edge_m(buffer(i-1):buffer(i))...
                                              + (n==1)*edge_delta(buffer(i-1):buffer(i))));
    exp.secondary_steps(i-1) = tau_seq_itr;
end

%% Save and exit
if exp.type == 0
    save('ExperimentMotor.mat','exp');
else
    save('ExperimentRudder.mat','exp');
end
