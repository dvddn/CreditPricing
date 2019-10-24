function price=KL_t(N,rho,recovery,ku,kd,I,k,nu)
% Function that computes the price of a mezzanine tranche thanks to the approximation of Kullback-Leibler
%
% INPUT: 
% N:        notional
% rho:      correlation between mortgages
% ku:       detachment point of the mezzanine tranche (upper notional bound)
% kd:       subordination point of the mezzanine tranche (lower notional bound)
% I:        numebr of mortgages
% k:        effective parameter: notional of each credit
% nu:       # of degrees of freedom of the t distribution
% 
% 
%OUTPUT:
% price:    price of a mezzanine tranche under the assumption of the KL approximation

u=ku/(1-recovery);
d=kd/(1-recovery);
L= @(z) min(max(z-d, 0), u-d)/(u-d);

% Default probability
p=@(eta) normcdf((eta./sqrt(1-rho)));

C1= @(z) sqrt(I./(2*pi.*(1-z).*z));
K= @(z, p) z.*log(z./p)+(1-z).*log((1-z)./(1-p));

% Mixing density
integrand = @(eta) integral(@(w) exp(-1./(2*rho).*(eta-k.*sqrt(w/nu)).^2).*w.^(nu/2-1).*exp(-w/2),0,100, 'ArrayValued', true);
density = @(eta) 1./(sqrt(rho*pi).*2^((nu+1)/2)*gamma(nu/2))*integrand(eta);   

g= @(z, eta) L(z).*density(eta).*exp(-I.*K(z, p(eta))).* C1(z);

loss= quadgk (@(eta) arrayfun( @(u) quadgk (@(z) g(z,u), 0, 1), eta), -15, 15);

% Price of the mezzanine tranche 
price=N*(ku-kd)*(1-loss);

end %function KL_t
