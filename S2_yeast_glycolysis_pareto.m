% Find equations for yeast glycolysis state variable 2

% define libarary parameters
laurentorder = 0;
polyorder = 6;
usesine = 0;
dyorder = 1;

% search for equation for the 2nd state-variable
% the pareto front and details for the 2nd variable were plotted in the paper
% clear results for other state variables
clear Theta Thetastring Xi indTheta lambdavec numterms errorv indopt
clear indTheta1 Xi1 numterms1 nT

% pool Data  (i.e., build library of nonlinear time series)
[Theta, Thetastring] = poolDatady(xt,n,polyorder,usesine, laurentorder, dxt(:,2), dyorder);

% compute Sparse regression using ADM
tol = 2e-3;
[Xi, indTheta, lambdavec, numterms, errorv] = ADMpareto(Theta, tol, plottag);

%find the optimal solution and print it out

indopt = find(min(abs(numterms(1:end-1).*errorv(1:end-1))) == abs(numterms(1:end-1).*errorv(1:end-1)));

Xi{indopt}
Thetastring(indTheta{indopt})'

save('2nd_state_variable_pareto.mat')

% Plot Pareto 
figure
semilogy(numterms, errorv, 'o')
xlabel('Number of terms')
ylabel('Validated RMSE')
title('Pareto Front for Validation')