function [para, ycdf] = fit_unchanged(x, y)


options = optimoptions('fmincon','MaxIterations',1e5,...
          'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);

xplot = 0:0.001:1.2;

% set up model parameters
mu = .45; % par
sigma = .05; % par
chance = 0.175; % para
asymt = 0.95; % constant
lb=[0 0 0 0]; ub = [1 inf 1 1];

alpha = 500; % regularization parameter
slope0 = .1; % slope prior

% input data
LL = @(params) -nansum(y.*log(params(3)+(params(4)-params(3))*normcdf(x,params(1),params(2)))...
        + (1-y).*log(1-(params(3)+(params(4)-params(3))*normcdf(x,params(1),params(2)))))...
        + alpha*(params(2)-slope0)^2; % control the
   % slope
    %  + alpha*sum(abs(params)); % lasso
    %+ alpha*(params*params'); % ridge
   % 
para = fmincon(LL,[mu sigma chance asymt],[],[],[],[],lb,ub,[],options);
ycdf = para(4)*normcdf(xplot,para(1),para(2))+para(3)*(1-normcdf(xplot,para(1),para(2)));