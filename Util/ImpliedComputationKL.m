function [rho, values, exitflag]= ImpliedComputationKL(Data,Market,model,nu)
% Function computing correlation and implied correlations for t student and
% dobule t student models for MBS pricing
%
% INPUTS:
% Data: Struct containing data about the MBS (detachment points,
%       notional,default probability)
% Market: Struct containing market data (implied correlations, prices)
% model: 0 for Double t model, 1 for Single t model
% nu: # of degrees of freedom of the model
%
% OUTPUTS:
% rho: correlation of the model
% values: if flag=1: Compound Implied correlation of the model
%         if flag=2: Base Implied correlation of the model
%         empty otherwise
% exitflag: In case of calibration error (rho=0 or implied correlation=0)
%           the exitflag is set to 0. Default value is 1.
%
% USES:
% KL_vasicek: computes the price of a MBS tranche through gaussian
%              vasicek model with KL approximation
% KL_double_t: computes the price of a MBS tranche through double student
%               t model with KL approximation
% KL_t: computes the price of a MBS tranche through student t model with
%        KL approximation

exitflag=1;
k2=tinv(Data.p,nu);
options = optimset('FunValCheck', 'on');
if model
    delta = @(x) KL_equity(Data.N,Market.base(1),Data.recovery,Data.ku(1),Data.I,Data.k1,2,nu)...
                    -KL_equity(Data.N,x,Data.recovery,Data.ku(1),Data.I,k2,1,nu);
else
    delta = @(x) KL_equity(Data.N,Market.base(1),Data.recovery,Data.ku(1),Data.I,Data.k1,2,nu)...
                    -KL_equity(Data.N,x,Data.recovery,Data.ku(1),Data.I,inversep(x,Data.p,nu),0,nu);
end

rho=fzero(delta, Market.base(1),options); %correlation of the model
k3=inversep(rho,Data.p,nu); %threshold k for double t model
values=zeros(4,1);
for j=1:4
        if model
            deltaimplied =@(x) KL_equity(Data.N,x,Data.recovery,Data.ku(j+1),Data.I,Data.k1,2,nu)...
                                -KL_equity(Data.N,rho,Data.recovery,Data.ku(j+1),Data.I,k2,1,nu);
        else
            deltaimplied =@(x) KL_equity(Data.N,x,Data.recovery,Data.ku(j+1),Data.I,Data.k1,2,nu)...
                                -KL_equity(Data.N,rho,Data.recovery,Data.ku(j+1),Data.I,k3,0,nu);
        end
    values(j)=fzero(deltaimplied, Market.base(j) ,options);
end

values=[Market.comp(1);values];

if ~(prod(values)&&rho)
    exitflag=0;
end

plotimplied(values,Market.base,Data.ku);



end %function ImpliedComputationKL