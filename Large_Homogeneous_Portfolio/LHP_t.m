function price=LHP_t(N,rho,recovery,ku,kd,k,nu)
% Price of a tranche with the closed formula for double t-Student model under HP
% assumptions
%
%INPUT
% N:           Notional of the portfolio
% rho:         Correlation
% recovery:    Recovery rate
% ku:          Detachment point of the tranche
% kd:          Attachment point of the tranche
% I:           Number of mortgages
% k:           Effective parameter (tinv(default probability of each mortgage,nu))
% nu:          Degrees of freedom of the t-Student distribution

u=ku/(1-recovery);
d=kd/(1-recovery);

% Loss function
loss_tranche=@(x) min(max(x-d,0),u-d)./(u-d);

integrand = @(x) integral(@(w) exp(-1./(2*rho).*(sqrt(1-rho)*norminv(x)-k.*sqrt(w/nu)).^2).*w.^(nu/2-1).*exp(-w/2),0,100,'ArrayValued', true);
density = @(x) sqrt(1-rho)./(normpdf(norminv(x)).*sqrt(rho*pi).*2^((nu+1)/2)*gamma(nu/2)).*integrand(x);
  
integrand2 = @(x) loss_tranche(x).*density(x);
loss=integral(@(x) arrayfun(integrand2, x) ,0,1, 'ArrayValued', true);

price=(ku-kd)*N*(1-loss);
 
end % function LHP_t