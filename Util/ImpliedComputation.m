function [rho, values, exitflag]= ImpliedComputation(Data,Market,model,flag,nu)
% Function computing correlation and implied correlations for t student and
% dobule t student models for MBS pricing
%
% INPUTS:
% Data: Struct containing data about the MBS (detachment points,
%       notional,default probability)
% Market: Struct containing market data (implied correlations, prices)
% flag: flag indicating the type of calibration to be used
%       1: Calibration through Compound Correlation 
%       2: Calibration through Base Correlation 
%       3: Calibration through Mezzanine Prices - for coherency checks
%       4: Calibration through Cumulated Equity Prices -  for coherency checks
% model: 0 for Double t model, 1 for Single t model
% nu: # of degrees of freedom of the model
%
% OUTPUTS:
% rho: correlation of the model
% values: if flag=1: Compound Implied correlation of the model
%         if flag=2: Base Implied correlation of the model
%         if flag=3: Mezzanine prices computed with the model
%         if flag=4: Cumulative equity prices computed with the model
% exitflag: In case of calibration error (rho=0 or implied correlation=0)
%           the exitflag is set to 0. Default value is 1.
%
% USES:
% LHP_vasicek: computes the price of a MBS tranche through gaussian
%              vasicek model with LHP assumption
% LHP_double_t: computes the price of a MBS tranche through double student
%               t model with LHP assumption
% LHP_t: computes the price of a MBS tranche through student t model with
%        LHP assumption

exitflag=1;
k2=tinv(Data.p,nu);
options = optimset('FunValCheck', 'on');

if model
    delta = @(x) LHP_vasicek(Data.N,Market.comp(1),Data.recovery,Data.ku(1),0,Data.k1)-LHP_t(Data.N,x,Data.recovery,Data.ku(1),0,k2,nu);
else
    delta = @(x) LHP_vasicek(Data.N,Market.comp(1),Data.recovery,Data.ku(1),0,Data.k1)...
                    -LHP_double_t(Data.N,x,Data.recovery,Data.ku(1),0,inversep(x,Data.p,nu),nu);
end

rho=fzero(delta, Market.base(1),options); %correlation of the model
k3=inversep(rho,Data.p,nu); % threshold k for double t model
values=zeros(4,1);

switch flag
    case 1
        for j=1:4
            if model
                deltaimplied =@(x) LHP_vasicek(Data.N,x,Data.recovery,Data.ku(j+1),Data.kd(j+1),Data.k1)...
                                    -LHP_t(Data.N,rho,Data.recovery,Data.ku(j+1),Data.kd(j+1),k2,nu);
            else
                deltaimplied =@(x) LHP_vasicek(Data.N,x,Data.recovery,Data.ku(j+1),Data.kd(j+1),Data.k1)...
                                    -LHP_double_t(Data.N,rho,Data.recovery,Data.ku(j+1),Data.kd(j+1),k3,nu);
            end
            values(j)=fzero(deltaimplied, Market.comp(j) ,options);
        end
        values=[Market.comp(1);values];
        if ~(prod(values)&&rho)
            exitflag=0;
        end
        plotimplied(values,Market.comp,Data.ku);
        
    case 2
        for j=1:4
            if model
                deltaimplied =@(x) LHP_vasicek(Data.N,x,Data.recovery,Data.ku(j+1),Data.kd(1),Data.k1)...
                                    -LHP_t(Data.N,rho,Data.recovery,Data.ku(j+1),Data.kd(1),k2,nu);
            else
                deltaimplied =@(x) LHP_vasicek(Data.N,x,Data.recovery,Data.ku(j+1),Data.kd(1),Data.k1)...
                                    -LHP_double_t(Data.N,rho,Data.recovery,Data.ku(j+1),Data.kd(1),k3,nu);
            end
            values(j)=fzero(deltaimplied, Market.base(j) ,options);
        end
        values=[Market.comp(1);values];
        if ~(prod(values)&&rho)
            exitflag=0;
        end
        plotimplied(values,Market.base,Data.ku);
        
        
    case 3
        for j=1:5
            if model
                values(j) = LHP_t(Data.N,rho,Data.recovery,Data.ku(j),Data.kd(j),k2,nu);
            else
                values(j) = LHP_double_t(Data.N,rho,Data.recovery,Data.ku(j),Data.kd(j),k3,nu);
            end
        end
        if ~rho
            exitflag=0;
        end
        plotimplied(values,Market.mezz,Data.ku);
        
        
    case 4
        for j=1:5
            if model
                values(j) = LHP_t(Data.N,rho,Data.recovery,Data.ku(j),Data.kd(1),k2,nu);
            else
                values(j) = LHP_double_t(Data.N,rho,Data.recovery,Data.ku(j),Data.kd(1),k3,nu);
            end
        end
        if ~rho
            exitflag=0;
        end
        plotimplied(values,Market.eq,Data.ku);
        
end

end