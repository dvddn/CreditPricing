function price=LHP_vasicek(N,rho,recovery,ku,kd,k)
% Price of a tranche with the closed formula for Vasicek model under LHP
% assumptions
%
% INPUTS:
% N:           Notional of the portfolio
% rho:         Correlation
% recovery:    Recovery rate
% ku:          Detachment point of the tranche
% kd:          Attachment point of the tranche
% k:           Effective parameter (norminv(default probability of each mortgage))

u=ku/(1-recovery);
d=kd/(1-recovery);

% Loss function
loss_tranche=@(x) min(max(x-d,0),u-d)./(u-d);

density=@(x) sqrt((1-rho)/rho)*exp(0.5*(norminv(x)).^2-(((k-norminv(x).*sqrt(1-rho)).^2)/(2*rho)));
loss=quadgk(@(x) loss_tranche(x).*density(x),0,1);

price=(ku-kd)*N*(1-loss);

end % function LHP_vasicek