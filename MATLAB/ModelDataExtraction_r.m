%load('ExperimentDataMotor.mat');
%load('ExperimentDataMotor.mat');
load('ActualExperimentDataRudder.mat');
%% EXTRACTION OF 'STEADY STATE'- AND DAMPING DATA
% Average steady state regions
x_ss_avg = zeros(2,length(tau_sequence));

for i = 1:length(x_ss_avg)
    x_ss_avg(:,i) = mean( (x(:,ss_timestamps(1,i):ss_timestamps(2,i))') );
end

% Extract damping measurements
sigma_r = tau_sequence(2,:);

%% EXTRACTION OF INERTIA DATA - SOG
% Declare result buffers for esimated inertias and transient response midrange values
inertia_r            = zeros(3,length(ss_timestamps)-1);
x_transient_midrange = zeros(2,length(ss_timestamps)-1);

% Declare buffers for simulation variables
r_sim        = NaN(1,length(x));
deltaTau_r_i = NaN(1,length(x));

% This for-loop iterates through all the transient responses and uses a
% binary search to identify an inertia value that best fits the response of
% a first order system to transient response in question.
max_attempts = 200;
error_quasi_convex = zeros(2,max_attempts);
for i = 1:(length(ss_timestamps) -1)
    transient_response_indices = ss_timestamps(2,i):ss_timestamps(1,i+1);
   
    % Set up simulation parameters for the i'th transient response
    r_sim(transient_response_indices(1)) = x(2,transient_response_indices(1)); % Initial SOG
    k_i = (sigma_r(i+1) - sigma_r(i))/(x_ss_avg(2,i+1) - x_ss_avg(2,i));   % Linearized damping term
    
    inertia_search_bounds = [0.001 100]; % We limit our search to be within these extrema
    error = -1;                        % Quadratic error term: (1/N)*sum((measured - estimated)^2)
    
    % Initialize logical variables
    attempt = 1;

    inertia_best_pos_fd = 0;   % Inertia tested with currently best error and finite difference > 0 
    inertia_best_neg_fd = 0;   % Inertia tested with currently best error and finite difference < 0 
    
    inertia_log_pos_fd = zeros(3,max_attempts)*Inf;    % Logging of inertia values,
    inertia_log_neg_fd = zeros(3,max_attempts)*(-Inf); % (1) m_U_i value, (2) m_U_i error
    
    % Search until error requirements are fullfilled, or max attempts reached.
    % We use a variaton of binary search to find a suitable value.
%     test_vals = linspace(0.01,50,500);
    while attempt < max_attempts %(error > 0.000001 || error < 0) && 
%         m_r_i = test_vals(attempt);
        error     = 0; % Error of simulation with  m_U_i
        error_hat = 0; % Error of simulation with 'm_U_i + delta'
        
        % Start at search bounds
        if attempt <= 2
            m_r_i = inertia_search_bounds(attempt);
        else
            m_r_i = (inertia_best_neg_fd + inertia_best_pos_fd)/2;
        end
        
        % This for-loop simulates the i'th transient response again, but
        % with the predicted inertia used in the last loop + a small
        % offset so that we can determine the sign of the error gradient
        for j = transient_response_indices
            deltaTau_r_i(j) = tau(2,j) - tau_sequence(2,i);                             
            r_sim(j+1) = r_sim(j) + h*(1/(m_r_i + 0.00001)*(deltaTau_r_i(j) + k_i*(x_ss_avg(2,i) - r_sim(j)))); % Euler integration of U, now with slightly different m_U_i
            error_hat = error_hat + (x(2,j)-r_sim(j))^2;
        end
        
        % This for-loop simulates the i'th transient response with a predicted inertia value
        for j = transient_response_indices
            deltaTau_r_i(j) = tau(2,j) - tau_sequence(2,i);                             
            r_sim(j+1) = r_sim(j) + h*(1/m_r_i*(deltaTau_r_i(j) + k_i*(x_ss_avg(2,i) - r_sim(j)))); % Euler integration of U
            error = error + (x(2,j)-r_sim(j))^2;
        end
        
        % Calculate errors, and evaluate and log inertia prediction 
        error     = 1/length(transient_response_indices) * error;
        error_hat = 1/length(transient_response_indices) * error_hat;
        
        error_fd  = error_hat - error; % Error finite difference, tells us something about the error gradient
        
        error_quasi_convex(1,attempt) = m_r_i; % For report
        error_quasi_convex(2,attempt) = error;
        
        % Working under the assumption that only small values of m_U_i will
        % result in error = NaN we can still make use of these values by
        if isnan(error)
            error = 1000 + m_r_i;   % We use '+ m_U_i' so that the newest estimates, allthough bad, are better than the old
            error_fd = -1;
        end
        
        if error_fd < 0
            inertia_log_neg_fd(:,attempt) = [m_r_i error error_fd]';
            inertia_log_neg_fd = sortrows(inertia_log_neg_fd',2)';
            
            inertia_best_neg_fd = m_r_i;
        else
            inertia_log_pos_fd(:,attempt) = [m_r_i error error_fd]';
            inertia_log_pos_fd = sortrows(inertia_log_pos_fd',2)';
            
            inertia_best_pos_fd = m_r_i;
        end
        
        attempt = attempt +1;
    end
    
    % Save results
    if abs(inertia_log_neg_fd(2,1)) < inertia_log_pos_fd(2,1)
        inertia_r(1,i) = inertia_log_neg_fd(1,1);
        inertia_r(2,i) = inertia_log_neg_fd(2,1);
        inertia_r(3,i) = attempt;   
    else
        inertia_r(1,i) = inertia_log_pos_fd(1,1);
        inertia_r(2,i) = inertia_log_pos_fd(2,1);
        inertia_r(3,i) = attempt;
    end
    x_transient_midrange(:,i) = [ (x_ss_avg(1,i+1) + x_ss_avg(1,i))/2;
                                  (x_ss_avg(2,i+1) + x_ss_avg(2,i))/2 ];
                              
%     figure(15);
%     clf(figure(15),'reset');
%     
%     [error_quasi_convex(1,:), indices] = sort(error_quasi_convex(1,:));
%     error_quasi_convex(2,:) = error_quasi_convex(2,indices);
%     hold on
%     grid on
%     plot(error_quasi_convex(1,:), error_quasi_convex(2,:),'.-r')
%     %ylim([0 0.01])
%     xlabel('$m_{r_i}$')
%     ylabel('$Error$')
end

%% POST PROCESSING
% Clean and mirror data
% - Damping data

sigma_r   = sigma_r;% remove_irrelevant_step_responses(sigma_r, motor_steps);
x_sigma_r = x_ss_avg;%[remove_irrelevant_step_responses(x_ss_avg(1,:), motor_steps);
             %remove_irrelevant_step_responses(x_ss_avg(2,:), motor_steps)];

sigma_r   = [sigma_r(:)'    -sigma_r(:)'];
x_sigma_r = [x_sigma_r(1,:)  x_sigma_r(1,:) ;
             x_sigma_r(2,:) -x_sigma_r(2,:)];

% - Inertia data
inertia_r   =  remove_irrelevant_step_responses(inertia_r(1,:), motor_steps);
x_inertia_r = [remove_irrelevant_step_responses(x_transient_midrange(1,:), motor_steps);
               remove_irrelevant_step_responses(x_transient_midrange(2,:), motor_steps)];

inertia_r   = [inertia_r(1,:)  inertia_r(1,:)];
x_inertia_r = [x_inertia_r(1,:)  x_inertia_r(1,:);
               x_inertia_r(2,:) -x_inertia_r(2,:)];

% Remove yaw rate offset
tau_sequence = [remove_irrelevant_step_responses(tau_sequence(1,:), motor_steps);
                remove_irrelevant_step_responses(tau_sequence(2,:), motor_steps)];
tau_sequence = [tau_sequence(1,:)  tau_sequence(1,:);
                tau_sequence(2,:) -tau_sequence(2,:)];
            
% for i = 1:length(tau_sequence)
%     if tau_sequence(2,i) == 0
%         x_sigma_r(2,i) = 0;
%         if i <= length(x_inertia_r) 
%             x_inertia_r(2,i) = 0;
%         end
%     end
% end


save('ActualMeasurementsRudder.mat','x_sigma_r','sigma_r','x_inertia_r', 'inertia_r');
%% PLOT
% Damping measurements
dampingMeasurements = figure(2);
clf(dampingMeasurements,'reset')

grid on
scatter3((x_sigma_r(2,:).*180/pi)',x_sigma_r(1,:)',sigma_r(:)','r')
set(gca,'Xdir','reverse')
xlabel('ROT [deg/s]')
ylabel('SOG [m/s]')
zlabel('$\sigma_r$')
xlim([-50 50])
ylim([0 20])
zlim([-1 1])
pbaspect([1 1 2])

%% Inertia measurements
inertiaSOGMeasurements = figure(8);
%clf(inertiaSOGMeasurements,'reset');

hold on
grid on

inertia_r_Tel = zeros(2,length(inertia_r));
for i = 1: length(inertia_r)
    
end

scatter3(x_inertia_r(1,:)',(x_inertia_r(2,:).*180/pi)',inertia_r(1,:)','.b')
xlabel('SOG [m/s]')
ylabel('ROT [deg/s]')
zlabel('$m_r$')
xlim([0 17])
ylim([-20 20])
zlim([0 20])
pbaspect([1 1 2])


%% Simulation vs real data
inertiaSOGsim = figure(4);
clf(inertiaSOGsim,'reset')

hold on
grid on

plot((1:length(r_sim)).*h, x(2,1:length(r_sim)).*180/pi,...
     (1:length(r_sim)).*h, r_sim.*180/pi)

xlabel('time [s]')
ylabel('ROT [deg/s]')

legend('r_{actual}', 'r_{approx}')
hold off

