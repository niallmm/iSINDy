% Michaelis menton kinetics

% Copyright 2017, All Rights Reserved
% Code by Niall Mangan for paper "Inferring biological networks by sparse
% identification of nonlinear dynamics"
% by N. M. Mangan S. L. Brunton, J. L. Proctor, and J. N. Kutz


clear all, close all, clc
figpath = '../figures/';
addpath('./utils');
addpath('./bioutils');

%% Define System

%size of system
n = 1;
% rate constants
Vmax = 1.5; % maximum reaction rate
Km = 0.3;   % half max reaction rate
jin = 0.6;  % influx of substrate

% reaction function
MMKinetics = @(x)jin-Vmax.*x./(Km+x);

%% generate Data
measure = 2; % number of initial conditions
dt = 0.1; % time step saved
tspan=[0:dt:4]; % time vector
N = length(tspan);

% intial condition of concentration
x0 = 0.5;

for ii = 1:measure
    % Integrate
    options = odeset('RelTol',1e-7,'AbsTol',1e-7);
    [t1,x1]=ode45(@(t,x)MMKinetics(x),tspan,x0,options);
    tt(:,ii) = t1;
    x(:,:,ii) = x1;
    
    x0 = 2*x0; % get a new initial condition
end
%% add noise & calculate derivative
eps = 1e-4; %magnitude of noise
xn = x + eps*randn(size(x)); % add normally distributed measuement error

figure(1)

%calculate exactly and add error
xt = [];dxt= []; t = [];

for ll =1:measure
    for ii=1:length(x)
        dxf(ii,:,ll) = MMKinetics(x(ii,:, ll));
    end
    epsdt =0;
    dxf = dxf+epsdt*randn(size(dxf));
    
    dxt = [dxt; dxf(:,:,ll)];
    xt = [xt; xn(:,:,ll)];
    t = [t; tt(:, ll)];
end

% % % use TVRegDiff to numerically calculate the derivative
% % % this is very sensative, and takes some playing around to make work
% % % may require more data or lower noise.
% mid = floor(length(tspan)/10);
% 
% xt =[]; t= []; dxt = [];
% clear dxf
% for mm = 1:measure %for each initial condition
%     xtemp = xn(:,:,mm);
%     
%     for nn =1:n % for each variable
%         
%         dx(:, nn) = TVRegDiff( xtemp(:,nn), 5, 1e4, [], 'small', 1e17, dt, 1, 1 );
%         dxf(:,nn, mm)= dx(:,nn); %store in tensor
%         % algorithm does better if you calculate the x values back from
%         % this numerically calculated derivative:
%         x2(:,nn) = cumsum(dx(:,nn))*dt;
%         x2(:,nn) = x2(:,nn) - (mean(x2(mid:end-mid,nn)) - mean(xtemp(mid:end-mid,nn)));
%         x2f(:,nn,mm) = x2(:,nn); % store in tensor 
%     end
%     % stack the time series values and derivatives into a single vector
%     % there are some issues with calculating a derivative numerically at
%     % the edges so not all points are kept
%     dxt = [dxt;dxf(4:end-1,:,mm)];
%     xt = [xt; x2f(4:end-1,:,mm)];
%     t = [t; tt(4:end, mm)];
% end


% Plot time series and derivatives of "training time series"
figure(5)
hold off
plot(t ,xt, 'o')
xlabel('time')
ylabel('concentrations')
title('training time series')

figure(6) 
hold off
plot(t, dxt, 'o')
xlabel('time')
ylabel('derivative of concentrations w/ time')
title('training derivative time series')

%% set the functions to search
% may not be the same as the functions we used to generate the data
laurentorder = 0; % n, for 1/x^n terms in the library
polyorder = 4; % n, for polynomial terms in the library x^n
usesine = 0; % adds sine to the library
dyorder = 1; % n for (dx/dt)^n terms in the library

% pool Data  (i.e., build library of nonlinear time series)
[Theta, Thetastring] = poolDatady(xt,n,polyorder,usesine, laurentorder, dxt, dyorder);

%% compute Sparse regression using ADM
pflag = 2; %plot output. 1 will plot some pareto front output, 2 will plot ADM algorithm details.
tol = 1e-5;
[Xi, indTheta, lambdavec, numterms, errorv] = ADMpareto(Theta, tol, pflag);

%% compare infered model and original model for new initial conditions

figure(12)
x0 = [0.1 0.6 1.2];
dt = 0.1;
tspan2=[0:dt:20];

for i = 1:length(x0)
    [t1,x1]=ode45(@(t,x)MMKinetics(x),tspan2,x0(i),options);
    x1 = x1+ eps*randn(size(x1));
    plot(t1,x1, 'o')
    hold on
end
%% numerically simulate time series for each SINDy discovered model to validate
[libsize, nummods] = size(Xi); % find the number of models

for kk = 1:nummods % loop through sparse coefficient vectors
    kk % display the set of coefficients we are on
    newdxdt = @(x)-Xi(1:size(Xi,1)/2, kk)'*[1 x x*x x*x*x x*x*x*x]'...
        /(Xi(size(Xi,1)/2+1:end, kk)'*[1 x x*x x*x*x x*x*x*x]');

    tspan=[0:dt:20];
    for i = 1:length(x0)
        [t2,x2]=ode23s(@(t,x)newdxdt(x), tspan, x0(i), options);
        hold on
        plot(t2,x2)
    end
    xlabel('time')
    ylabel('x')
    drawnow
    hold off
    
    % save validation time series
    x2val{kk} = x2;
    t2val{kk} = t2;
    
    % store Root Mean Square error between data and validation time series
    RMSE(kk) = sqrt(mean((x2-x1(length(x2),1)).^2));
    
end

% Plot Pareto with validation step
figure
semilogy(numterms, RMSE, 'o')
xlabel('Number of terms')
ylabel('Validated RMSE')
title('Pareto Front for Validation')
