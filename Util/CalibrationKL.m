function [nu] = CalibrationKL(Data,Market,model,tol,maxIter,int)
% function that calibrates either Single t or Double t models for MBS
% pricing with through implied correlations (compound or base) or prices
% (cumulatives or mezzanines) through KL approx
% 
% INPUTS:
% Data:     Struct containing data about the MBS (detachment points, notional,default probability)
% Market:   Struct containing market data (implied correlations, prices)
% model:    0 for Double t model, 1 for Single t model
% tol:      error tolerance of the calibration algorithm
% maxIter:  max number of iterations of the calibration algorithm
% int:      search interval for nu, containing two numbers >2 in increasing order
% 
% OUTPUTS:
% nu:       # of degrees of freedom of the model
% rho:      correlation of the model
% values:   if flag=1: Compound Implied correlation of the model
%           if flag=2: Base Implied correlation of the model. 
%           Empty otherwise
% exitflag: In case of calibration error (rho=0 or implied correlation=0) the exitflag is set to 0. Default value is 1.
% 
% USES:
% KL_vasicek:  computes the price of a MBS tranche through gaussian vasicek model with KL approximation
% KL_double_t: computes the price of a MBS tranche through double student t model with KL approximation
% KL_t:        computes the price of a MBS tranche through student t model with KL approximation
% MSEevalKL:  computes Mean Square Error of the implied base correlation for single or double t models

if nargin<5
    tol=1e-4;  maxIter=100;  int=[2+5*model,20-8*model];
elseif nargin<6
    maxIter=100;  int=[2+5*model,20-8*model];
elseif nargin<7
    int=[2+5*model,20-8*model];
end

if (int(1)<2 || int(2)<int(1))
    error('Invalid Search Interval. Interval should be increasing and >2')
end


f=@(x) MSEevalKL(Data,Market.base,x,model);

options=optimset('Display', 'iter', 'MaxIter', maxIter,'TolX', tol);
nu=fminbnd(f,int(1),int(2),options);

end %function CalibrationKL