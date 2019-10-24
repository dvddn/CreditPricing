function [MSE]= MSEpricecomp(Data,CompoundCorr,Mezzanines,nu,flag)
% computes Mean Square Error of the mezzanine tranche prices for single or double t models
% 
% INPUTS:
% Data:         Struct containing data about the MBS (detachment points, notional,default probability)
% CompoundCorr: Market implied compound correlation for the 0-3% tranche
% Equity:       Market mezzanine tranche prices
% nu:           # of degrees of freedom of the model
% flag:         =0 for double t model, =1 for single t model
% 
% OUTPUTS:
% MSE:          Mean Square Error of the mezzanine tranche prices
% 
% USES:
% LHP_vasicek:  computes the price of a MBS tranche through gaussian vasicek model with LHP assumption
% LHP_double_t: computes the price of a MBS tranche through double student t model with LHP assumption
% LHP_t:        computes the price of a MBS tranche through student t model with LHP assumption


deltaimplied=zeros(4,1);    
k2=tinv(Data.p,nu);

if flag
    delta2 = @(x) LHP_vasicek(Data.N,CompoundCorr,Data.recovery,Data.ku(1),Data.kd(1),Data.k1) ...
                    -LHP_t(Data.N,x,Data.recovery,Data.ku(1),Data.kd(1),k2,nu);
else
    delta2 = @(x) LHP_vasicek(Data.N,CompoundCorr,Data.recovery,Data.ku(1),Data.kd(1),Data.k1) ...
                    -LHP_double_t(Data.N,x,Data.recovery,Data.ku(1),Data.kd(1),inversep(x,Data.p,nu),nu);
end

options = optimset('FunValCheck', 'on');

    try
        rho=fzero(delta2, 0.2,options);
        k3=inversep(rho,Data.p,nu);
        if isnan(rho)
            rho=0;
        end
    catch
        rho=0;
    end

    for j=1:4
        
        if flag
            deltaimplied(j) = Mezzanines(j+1)...
                                -LHP_t(Data.N,rho,Data.recovery,Data.ku(j+1),Data.kd(j+1),k2,nu);
        else
            deltaimplied(j) = Mezzanines(j+1)...
                                -LHP_double_t(Data.N,rho,Data.recovery,Data.ku(j+1),Data.kd(j+1),k3,nu);
        end
      
    end

MSE=sum(deltaimplied.^2)/5;

end %function MSEpricecomp