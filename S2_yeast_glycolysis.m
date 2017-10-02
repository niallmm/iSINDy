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

lambda = 2e-3;

jj = 1; % counter
num= 1; % initialize the number of nonzero terms found for the lambda
errorvec= 0;

MaxIter = 1e3;

% for now calculate null space using null function
nT = null(Theta);

[indTheta1, Xi1, numterms1] = ADMinitvary(nT,lambda,MaxIter,tol, plottag);

Thetastring(indTheta1)'
n0Xi = Xi1(Xi1~=0); % terms need to be rearranged to recover coefficients
n0Xi/n0Xi(end)

save('Results/2nd_state_variable_1.mat')
