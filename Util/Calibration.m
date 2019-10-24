function [nu] = Calibration(Data,Market,flag,model,tol,maxIter,int)
% function that calibrates either Single t or Double t models for MBS
% pricing with through implied correlations (compound or base) or prices
% (cumulatives or mezzanines)
% 
% INPUTS:
% Data:     Struct containing data about the MBS (detachment points, notional,default probability)
% Market:   Struct containing market data (implied correlations, prices)
% flag:     flag indicating the type of calibration to be used
%           1: Calibration through Compound Correlation
%           2: Calibration through Base Correlation
%           3: Calibration through Mezzanine Prices - just for double checking
%           4: Calibration through Cumulated Equity Prices - just for double checking
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
% LHP_vasicek:  computes the price of a MBS tranche through gaussian vasicek model with LHP assumption
% LHP_double_t: computes the price of a MBS tranche through double student t model with LHP assumption
% LHP_t:        computes the price of a MBS tranche through student t model with LHP assumption
% MSEevalcomp:  computes Mean Square Error of the implied compound correlation for single or double t models
% MSEevalbase:  computes Mean Square Error of the implied base correlation for the single or double t models
% MSEpricecomp: computes the Mean Square Error of the mezzanine prices for single or double t models
% MSEpricebase: computes the Mean Square Error of the cumulative prices for the single or double t models
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


switch flag
    case 1
        f=@(x) MSEevalcomp(Data,Market.comp,x,model);
    case 2
        f=@(x) MSEevalbase(Data,Market.base,x,model);
    case 3
        f=@(x) MSEpricecomp(Data,Market.base(1),Market.mezz,x,model);
    case 4 
        f=@(x) MSEpricebase(Data,Market.base(1),Market.eq,x,model);
end

options=optimset('Display', 'iter', 'MaxIter', maxIter,'TolX', tol);
nu=fminbnd(f,int(1),int(2),options);

end %function Calibration