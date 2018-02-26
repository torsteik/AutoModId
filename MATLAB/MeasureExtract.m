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
sigma = exp.tau_seq(n,:);

% Get steady state averages
x_ss = zeros(2,length(sigma));

for i = 1:length(sigma)
    if ss_timestamps(1,i) == 1
        x_ss(:,i) = exp.x(:,ss_timestamps(2,i));
    else
        x_ss(:,i) = mean((exp.x(:,ss_timestamps(1,i):ss_timestamps(2,i))'));
    end
end
%% Inertia measurements
m    = zeros(1,length(ss_timestamps)-1); % Inertia
x_ts = zeros(2,length(ss_timestamps)-1);

plot_bool   = input('\nTurn on intermediate plots?\n1 - Yes\n2 - No\nAnswer: ');
first_trans = input('\nSkip to specific transient?\nTransient number(1 = No): ');

% Go through all identified transients
for trans_itr = first_trans:length(ss_timestamps)-1
    % Skip invalid data and steps in secondary input
    if (ss_timestamps(1,trans_itr+1) == 1) || any(ismember(exp.secondary_steps,trans_itr))
        continue;
    end
    
    %% Set up
    % Create struct to store info about the transient
    trans     = struct;
    trans.id  = trans_itr;
    trans.n   = n;
    trans.h   = exp.h;
    trans.t   = ss_timestamps(2,trans.id):ss_timestamps(1,trans.id+1);
    trans.x   = exp.x(n,trans.t);
    trans.tau = exp.tau(n,trans.t);
    trans.ss  = [x_ss(n,trans.id);
                 x_ss(n,trans.id+1)];
    trans.tau_prev = exp.tau_seq(n,trans.id);
    trans.k   = (sigma(trans.id+1) - sigma(trans.id))/(x_ss(n,trans.id+1) - x_ss(n,trans.id)); % Linearized damping term
    
    trans.plot_bool = plot_bool;
    if trans.plot_bool == 1
        %% Plot
        r2d = 180/pi;
        
        trans.ts_fig = figure(301);
        clf(trans.ts_fig,'reset')
        set(trans.ts_fig,'defaultAxesColorOrder',[[0 0 1]; [1 0 0]]);
        
        trans.ts_ax = axes;
        xlabel(trans.ts_ax,'time [s]');
        hold(trans.ts_ax,'on');
        trans.ts_ax.XGrid = 'on';
        trans.ts_ax.YGrid = 'on';
        trans.ts_ax.LineWidth = 1;
        
        trans.ts_plot = plot(trans.ts_ax, trans.t.*trans.h, trans.x.*(r2d^(n-1)));
        trans.ts_ax.XLim = [trans.t(1) trans.t(end)].*trans.h;
        
        input('Press ENTER to continue')
    end
    
    %% Find an inertia measurement
    [m(trans.id),bfgs] = BFGS(trans);    
    x_ts(:,trans.id) = [ (x_ss(1,trans.id+1) + x_ss(1,trans.id))/2;
                         (x_ss(2,trans.id+1) + x_ss(2,trans.id))/2]; 

    if trans.plot_bool == 1
        %% Plot
        [~,x_sim] = eval_err(m(trans.id),trans);
        trans.ts_plot = plot(trans.ts_ax, trans.t.*trans.h, x_sim*(r2d^(trans.n-1)));
        input('Press ENTER to continue')
    end
end

%% Post processing
% Remove the invalid elements
invalid_ss = ismember(ss_timestamps(1,:),1);
invalid_ts = circshift(ismember(ss_timestamps(1,:),1),-1);

invalid_ts(exp.secondary_steps) = true;
invalid_ts(end) = true;

sigma(invalid_ss)  = [];
x_ss(:,invalid_ss) = [];

m(invalid_ts(1:end-1))      = [];
x_ts(:,invalid_ts(1:end-1)) = [];

% Assume symmetric rudder response
sigma = [sigma      sigma.*((-1)^(n-1))];
x_ss  = [x_ss(1,:)  x_ss(1,:);
         x_ss(2,:) -x_ss(2,:)];
     
m    = [m          m];
x_ts = [x_ts(1,:)  x_ts(1,:);
        x_ts(2,:) -x_ts(2,:)];
         
tau_seq = [exp.tau_seq(1,:)  exp.tau_seq(1,:);
           exp.tau_seq(2,:) -exp.tau_seq(2,:)];

invalid_ss = [invalid_ss invalid_ss];
invalid_ts = [invalid_ts invalid_ts];
%% Plot
set(groot, 'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

r2d = 180/pi;
% Plot #1: Damping/Sigma
sig_fig = figure(302+n);
clf(sig_fig,'reset')

sig_ax = axes;
sig_ax.PlotBoxAspectRatio = [1 1 1];
hold(sig_ax,'on');
sig_ax.XGrid = 'on';
sig_ax.YGrid = 'on';
sig_ax.ZGrid = 'on';
sig_ax.LineWidth = 1;


sig_scat = scatter3(sig_ax,x_ss(2,:).*r2d,x_ss(1,:),sigma);
sig_scat.Marker = '.';
sig_scat.MarkerEdgeColor = 'b';
sig_scat.MarkerEdgeAlpha = 0.9;
sig_scat.SizeData = 18;

xlabel(sig_ax,'ROT [deg/s]');
ylabel(sig_ax,'SOG [m/s]');
zlabel(sig_ax,strcat('$\sigma',(n==1)*'_U$',(n==2)*'_r$'));

sig_ax.XLim = [min(x_ss(2,:)) max(x_ss(2,:))].*r2d;
sig_ax.YLim = [0 max(x_ss(1,:))];
sig_ax.ZLim = [min(sigma) max(sigma)];

% Plot #2: Inertia/m
m_fig = figure(304+n);
clf(m_fig,'reset')

m_ax = axes;
m_ax.PlotBoxAspectRatio = [1 1 1];
hold(m_ax,'on');
m_ax.XGrid = 'on';
m_ax.YGrid = 'on';
m_ax.ZGrid = 'on';
m_ax.LineWidth = 1;


sig_scat = scatter3(m_ax,x_ts(2,:).*r2d,x_ts(1,:),m);
sig_scat.Marker = '.';
sig_scat.MarkerEdgeColor = 'b';
sig_scat.MarkerEdgeAlpha = 0.9;
sig_scat.SizeData = 18;

xlabel(m_ax,'ROT [deg/s]');
ylabel(m_ax,'SOG [m/s]');
zlabel(m_ax,strcat('$m',(n==1)*'_U$',(n==2)*'_r$'));

m_ax.XLim = [min(x_ts(2,:)) max(x_ts(2,:))].*r2d;
m_ax.YLim = [0 max(x_ts(1,:))];
m_ax.ZLim = [min(m) max(m)];

%% Save and exit
if n == 1
    save('ActualMeasurementsMotor.mat','sigma','m','x_ss','x_ts','invalid_ss','invalid_ts','tau_seq');
else
    save('ActualMeasurementsRudder.mat','sigma','m','x_ss','x_ts','invalid_ss','invalid_ts','tau_seq');
end

%% Remove data with tau_m > 0.6
