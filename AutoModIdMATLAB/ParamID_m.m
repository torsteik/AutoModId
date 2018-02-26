% SCRIPT #5
n = input('\nDo you want to load:\n1 - Motor measurements\n2 - Rudder measurements\nAnswer: ');
if n == 1
    load('ActualMeasurementsMotor.mat');
elseif n == 2
    load('ActualMeasurementsRudder.mat');
else
    'Invalid input'
end

plot_bool = input('\nTurn on intermediate plots?\n1 - Yes\n2 - No\nAnswer: ');
%% Preprocess measurements
% Normalize state measurements
[x_ss_norm, x_ss_mean, x_ss_stdev] = zscore(x_ss');
[x_ts_norm, x_ts_mean, x_ts_stdev] = zscore(x_ts');

x_ss_norm = x_ss_norm';
x_ts_norm = x_ts_norm';

%% Regression setup
% A linear regression is used for the identification.
% The 'linear regression with L1 regularization'-solver solves:
%   minimize ||phi*beta-y||^2 + lambda*sum|beta_i|
%
% It is called as follows:
%   [beta,status,history] = l1_ls(phi,y,lambda [,tar_gap[,quiet]])

x = x_ts_norm;
y = m;

% Basis function:
% - The basis function, phi(x), in the regression of m(x) is a
%   4th order polynomial. 
phi = zeros(16,length(x));

for i = 1:size(phi,2)
    phi(:,i) = [ 1         ... 
                 x(1,i)^1  0     *    x(2,i)^1  ...
                 x(1,i)^2  0*x(1,i)^1*x(2,i)^1           x(2,i)^2  ...
                 x(1,i)^3  0*x(1,i)^2*x(2,i)^1  x(1,i)^1*x(2,i)^2  0     *    x(2,i)^3  ...
                 x(1,i)^4  0*x(1,i)^3*x(2,i)^1  x(1,i)^2*x(2,i)^2  0*x(1,i)^1*x(2,i)^3  x(2,i)^4 ...
                 0];% The last term will be found using cross validation
end



thing1 = ~invalid_ts & tau_seq(n,:) | circshift(~invalid_ts,1);
tau_seq2 = tau_seq(:,thing1);
inputs = sort(unique(tau_seq2(n,:)));
ind1 = inputs(abs(inputs) ~= max(abs(inputs)));
ind2 = inputs(abs(inputs) ~= min(abs(inputs)));
ind3 = unique(tau_seq(3-n,:));

tau_seq3 = tau_seq(:,~invalid_ts);

hash = zeros(length(ind1),length(ind2),length(ind3));
j = 0;
for i = 1:length(tau_seq3(n,:))-1
    j = j + 1;
    while invalid_ts(j) == 1
        j = j + 1;
    end
    tau_start = tau_seq(n,j);
    tau_end   = tau_seq(n,j+1);
    tau_secondary = tau_seq(3-n,j);
    
    if sign(tau_start) < 0
        signum = -1;
    else
        signum = 1;
    end
    
    hash1 = find(ind1 == signum*min(abs(tau_start),abs(tau_end)));
    hash2 = find(ind2 == signum*max(abs(tau_start),abs(tau_end)));
    hash3 = find(ind3 == tau_secondary);
    
    hash(hash1,hash2,hash3) = 1;
end

K      = length(find(hash));                         % K-fold cv
folds  = crossvalind('Kfold',length(find(hash)), K);
hash(find(hash)) = folds; % Bytt ut med hash == 1 eler noe s�nt

folds = zeros(1,length(x_ts));
j = 0;
for i = 1:length(tau_seq3)
    j = j + 1;
    while invalid_ts(j) == 1
        j = j + 1;
    end
    tau_start = tau_seq(n,j);
    tau_end   = tau_seq(n,j+1);
    tau_secondary = tau_seq(3-n,j);
    
    if sign(tau_start) < 0
        signum = -1;
    else
        signum = 1;
    end
    hash1 = find(ind1 == signum*min(abs(tau_start),abs(tau_end)));
    hash2 = find(ind2 == signum*max(abs(tau_start),abs(tau_end)));
    hash3 = find(ind3 == tau_secondary);
    
    folds(i) = hash(hash1,hash2,hash3);
end

% Weighting scheme
weights = ones(1,length(x));

% for i = 1:length(x)
%     pair = find(folds == folds(i));
%     if length(pair) == 2
%         weights(i) = 1/(abs(m(pair(1)) - m(pair(2))))^2;
%     else
%         weights(i) = 0;
%     end
% end

for i = 1:length(x)
    pair = find(folds == folds(i));
    if length(pair) == 2
        weights(i) = 1/((m(pair(1)) + m(pair(2)))/2)^2;
    else
        weights(i) = 0;
    end
end

figure(1)
hold on
grid on
scatter3(180/pi*(x(2,:).*x_ts_stdev(2) + x_ts_mean(2)), (x(1,:).*x_ts_stdev(1) + x_ts_mean(1)),weights, 12, weights,'filled')
xlabel('ROT')
ylabel('SOG')
zlabel('Weight')


weights = diag(weights);

% Group measurements generated by equal inputs, tau, in cv
% (Dear Lord, have mercy on my soul.)
% - Create a hash table to tell which transients are valid           
thing1 = ~invalid_ts & tau_seq(n,:) | circshift(~invalid_ts,1);
tau_seq2 = tau_seq(:,thing1);
inputs = sort(unique(tau_seq2(n,:)));
ind1 = inputs(abs(inputs) ~= max(abs(inputs)));
ind2 = inputs(abs(inputs) ~= min(abs(inputs)));
ind3 = unique(tau_seq(3-n,:));

tau_seq3 = tau_seq(:,~invalid_ts);

hash = zeros(length(ind1),length(ind2),length(ind3));
j = 0;
for i = 1:length(tau_seq3(n,:))
    j = j + 1;
    while invalid_ts(j) == 1
        j = j + 1;
    end
    tau_start = tau_seq(n,j);
    tau_end   = tau_seq(n,j+1);
    tau_secondary = tau_seq(3-n,j);
    
    if sign(tau_start) < 0
        signum = -1;
    else
        signum = 1;
    end
    
    hash1 = find(ind1 == signum*min(abs(tau_start),abs(tau_end)));
    hash2 = find(ind2 == signum*max(abs(tau_start),abs(tau_end)));
    hash3 = find(ind3 == tau_secondary);
    
    hash(hash1,hash2,hash3) = 1;
end

K      = 10;                         % K-fold cv
folds  = crossvalind('Kfold',length(find(hash)), K);
hash(find(hash)) = folds; % Bytt ut med hash == 1 eler noe s�nt

folds = zeros(1,length(x_ts));
j = 0;
for i = 1:length(tau_seq3)
    j = j + 1;
    while invalid_ts(j) == 1
        j = j + 1;
    end
    tau_start = tau_seq(n,j);
    tau_end   = tau_seq(n,j+1);
    tau_secondary = tau_seq(3-n,j);
    
    if sign(tau_start) < 0
        signum = -1;
    else
        signum = 1;
    end
    hash1 = find(ind1 == signum*min(abs(tau_start),abs(tau_end)));
    hash2 = find(ind2 == signum*max(abs(tau_start),abs(tau_end)));
    hash3 = find(ind3 == tau_secondary);
    
    folds(i) = hash(hash1,hash2,hash3);
end

cv_ind = folds;
figure(2)
hold on
grid on
for i = 1:K 
    scatter3(x_ts(2,cv_ind == i), x_ts(1,cv_ind == i), m(cv_ind == i));
end


%% Hyperparameter #1 and #2: Assymptotic parameters,a and b
a = linspace(10,35,100/4);
b = linspace(-1.5,-1,30/2);

lambda_ab   = 5e-1;
err_ab = zeros(length(a),length(b));
for i = 1:length(a)
    for j = 1:length(b)
        % Fill in the basis function element using the current set of hyperparameters
        phi(end,:) = tanh(a(i)*(x(1,:) - b(j)));
        
        for fold = 1:K
            % Determine which measurements to use in validation in this iteration
            y_valid = y(folds == fold);
            y_train = y(folds ~= fold);

            phi_valid = phi(:, folds == fold);
            phi_train = phi(:, folds ~= fold);
            
            weights_train = weights(folds ~= fold, folds ~= fold);

            % Solve for the parameter vector, beta, using training set and current lambda
            beta = l1_ls(sqrt(weights_train)*phi_train',sqrt(weights_train)*y_train',lambda_ab,[],false);

            % Evaluate validation set with the obtained beta and store validation error
            err_ab(i,j) = err_ab(i,j) + sum( (beta'*phi_valid - y_valid).^2 )/length(y_valid);
        end
    end
end

cvErrors_ab = figure(18);
clf(cvErrors_ab,'reset');

[X, Y] = meshgrid(a,b);

hold on
grid on
surf(X,Y,err_ab')
hold off

% Choose the a,b - pair that gave the lowest cross validation error and
% calculate the corresponding values for the asymptotic terms in phi_inertia
[min_error_per_column,row_index_per_column] = min(err_ab);
[~,column_index]                            = min(min_error_per_column);

a = a(row_index_per_column(column_index));
b = b(column_index);

phi(end,:) = tanh(a*(x(1,:) - b));

%% Hyperparameter #3: regularization coefficient, lambda
% Perform CV
lambda_res = 100;                    % 'Resolution' of the cv
lambdas = linspace(0.01,2,lambda_res);

err_lambda = zeros(1,length(lambdas));
for i = 1:length(lambdas)
    for fold = 1:K
        % Use measurments belonging to fold for validation, train on rest
        y_valid = y(folds == fold);
        y_train = y(folds ~= fold);

        phi_valid = phi(:, folds == fold);
        phi_train = phi(:, folds ~= fold);
        
        weights_train = weights(folds ~= fold, folds ~= fold);
        
        % Solve for the parameter vector, beta, using training set and current lambda
        beta = l1_ls(sqrt(weights_train)*phi_train', sqrt(weights_train)*y_train', lambdas(i), [],true);
        
        % Evaluate validation set with the obtained beta and store validation error
        err_lambda(i) = err_lambda(i) + sum( (beta'*phi_valid - y_valid).^2 )/length(y_valid);
    end
end

% Choose the lambda value that gave the lowest validation error
[~,ind] = min(err_lambda);
lambda  = lambdas(ind);

if plot_bool == 1
    %% Plot
    err_lam_fig = figure(404+n);
    clf(err_lam_fig,'reset')

    err_lam_ax = axes;
    hold(err_lam_ax,'on');
    err_lam_ax.XGrid = 'on';
    err_lam_ax.YGrid = 'on';
    err_lam_ax.ZGrid = 'on';
    err_lam_ax.LineWidth = 1;

    err_lam_plot = plot(err_lam_ax, lambdas, err_lambda);
    err_lam_plot.LineWidth   = 1;
    
    xlabel(err_lam_ax,'$\lambda$');
    ylabel(err_lam_ax,'Error');
end

%% Complete the identification using the chosen lambda
beta = l1_ls(sqrt(weights)*phi',sqrt(weights)*y',lambda,[],true);

% Store results
model_m = struct;
model_m.phi  = phi;
model_m.beta = beta;
model_m.x_mean  = x_ts_mean;
model_m.x_stdev = x_ts_stdev;


%% Plot
r2d = 180/pi;
% Set up for surface plot
xy_surf = [linspace(min(x(1,:)),max(x(1,:)),30);
           linspace(min(x(2,:)),max(x(2,:)),30)];
z_surf  = zeros(length(xy_surf(1,:)),length(xy_surf(2,:)));
      
[X,Y] = meshgrid(xy_surf(1,:),xy_surf(2,:));
phi   = zeros(16,length(xy_surf(1,:)),length(xy_surf(2,:)));

for i = 1:length(xy_surf(1,:))
    for j = 1:length(xy_surf(1,:))
        phi(:,i,j) = [ 1         ... 
                       X(i,j)^1  0     *    Y(i,j)^1  ...
                       X(i,j)^2  0*X(i,j)^1*Y(i,j)^1          Y(i,j)^2  ...
                       X(i,j)^3  0*X(i,j)^2*Y(i,j)^1  X(i,j)^1*Y(i,j)^2  0     *    Y(i,j)^3  ...
                       X(i,j)^4  0*X(i,j)^3*Y(i,j)^1  X(i,j)^2*Y(i,j)^2  0*X(i,j)^1*Y(i,j)^3  Y(i,j)^4 ...
                       tanh(a*(X(i,j) - b))];
    end
end

for i = 1:length(xy_surf(1,:))
    for j = 1:length(xy_surf(2,:))
        z_surf(i,j) = phi(:,i,j)'*model_m.beta;
    end
end

X_unorm = model_m.x_mean(1) + model_m.x_stdev(1).*X;
Y_unorm =(model_m.x_mean(2) + model_m.x_stdev(2).*Y).*r2d;

% Plot modelsurface against measurements
set(groot, 'defaulttextinterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

m_id_fig = figure(406+n);
clf(m_id_fig,'reset')

m_id_ax = axes;
m_id_ax.PlotBoxAspectRatio = [1 1 1];
hold(m_id_ax,'on');
m_id_ax.XGrid = 'on';
m_id_ax.YGrid = 'on';
m_id_ax.ZGrid = 'on';
m_id_ax.LineWidth = 1;

m_id_scat = scatter3(m_id_ax,x_ts(2,:).*r2d,x_ts(1,:),y);
m_id_scat.Marker = 'o';
m_id_scat.MarkerEdgeColor = 'r';
m_id_scat.MarkerFaceColor = 'r';
m_id_scat.SizeData = 7;

m_id_surf = surf(m_id_ax, Y_unorm, X_unorm, z_surf);
m_id_surf.EdgeColor = 'interp';
m_id_surf.FaceColor = [0.6 0.6 0.6];
m_id_surf.LineWidth = 1;

xlabel(m_id_ax,'ROT [deg/s]');
ylabel(m_id_ax,'SOG [m/s]');
zlabel(m_id_ax,strcat('$m',(n==1)*'_U$',(n==2)*'_r$'));

m_leg = legend(m_id_ax, '$\mathbf{D}_{m}$');

%sig_id_ax.XLim = [min(x_ts(2,:)) max(x_ts(2,:))].*r2d;
%sig_id_ax.YLim = [0 max(x_ts(1,:))];
m_id_ax.ZLim = [min(y) max(y)];
