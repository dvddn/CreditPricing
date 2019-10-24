function price=KL_double_t(N,rho,recovery,ku,kd,I,k,nu)
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
% OUTPUT:
% price:    price of a mezzanine tranche under the assumption of the KL approximation

u=ku/(1-recovery);
d=kd/(1-recovery);
L= @(z) min(max(z-d, 0), u-d)/(u-d);

py= @(y) tcdf((k-sqrt(rho*(nu-2)/nu).*y)/(sqrt((1-rho)*(nu-2)/nu)),nu);

C1= @(z) sqrt(I./(2*pi.*(1-z).*z));
K= @(z, p) z.*log(z./p)+(1-z).*log((1-z)./(1-p));
g= @(z, y) L(z).*tpdf(y,nu).*exp(-I.*K(z, py(y))).* C1(z);

loss= quadgk (@(y) arrayfun( @(u) quadgk (@(z) g(z,u), 0, 1), y), -15, 15);
price=(1-loss)*N*(ku-kd);


end