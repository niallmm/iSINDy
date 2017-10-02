% Find equations for yeast glycolysis state variable 3

% define libarary parameters
laurentorder = 0;
polyorder = 3;
usesine = 0;
dyorder = 1;

% clear variables from other solutions.
clear Theta Thetastring Xi indTheta lambdavec numterms errorv indopt
clear indTheta1 Xi1 numterms1 nT

% pool Data  (i.e., build library of nonlinear time series)
[Theta, Thetastring] = poolDatady(xt,n,polyorder,usesine, laurentorder, dxt(:,3), dyorder);
% %initial lambda value, which is the value used for soft thresholding in ADM

tol = 2e-3;
lambda = 8e-3;

jj = 1; % counter
num= 1; % initialize the number of nonzero terms found for the lambda
errorvec= 0;

MaxIter = 1e3;

% for now calculate null space using null function
nT = null(Theta);

[indTheta1, Xi1, numterms1] = ADMinitvary(nT,lambda,MaxIter,tol, plottag);
Thetastring(indTheta1)'
n0Xi = Xi1(Xi1~=0); 
n0Xi/n0Xi(end)
save('Results/3rd_state_variable.mat')
