function price=LHP_double_t(N,rho,recovery,ku,kd,k,nu)
% Price of a tranche with the closed formula for double t-Student model under LHP
% assumptions
%
%INPUT
% N:           Notional of the portfolio
% rho:         Correlation
% recovery:    Recovery rate
% ku:          Detachment point of the tranche
% kd:          Attachment point of the tranche
% k:           Effective parameter (tinv(default probability of each mortgage,nu))
% nu:          Degrees of freedom of the t-Student distribution

u=ku/(1-recovery);
d=kd/(1-recovery);

% Loss function
loss_tranche=@(x) min(max(x-d,0),u-d)./(u-d);

density=@(x) sqrt((1-rho)/rho)*tpdf((tinv(x,nu)*sqrt(1-rho)*sqrt((nu-2)/nu)-k)/sqrt(rho*(nu-2)/nu),nu)./(tpdf(tinv(x,nu),nu));
loss=quadgk(@(x) loss_tranche(x).*density(x),0,1);

price=(ku-kd)*N*(1-loss);

end % function LHP_double_t