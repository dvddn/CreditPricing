function price=KL_equity(N,rho,recovery,ku,I,k,model,nu)
% Function that computes the price of a mezzanine tranche thanks to the approximation of Kullback-Leibler
%
% INPUT: 
% N:        notional
% rho:      correlation between mortgages
% ku:       detachment point of the mezzanine tranche (upper notional bound)
% I:        numebr of mortgages
% k:        effective parameter: defaul threshold of each credit
% model:    =0 for double t
%           =1 for single t
%           =2 for vasicek
% nu:       # of degrees of freedom of the t distribution
% 
% OUTPUT:
% price:    price of a mezzanine tranche under the assumption of the KL approximation
% 
% USES:
% KL_double_t:  computes price of a mezzanine through KL approx with double t model
% KL_t:         computes price of a mezzanine through KL approx with student t model         
% KL_vasicek:   computes price of a mezzanine through KL approx with vasicek


if (nargin==7 && model==2)
    nu=0;
end
    
switch model
    case 0
        mezzanine=N*(1-ku)-KL_double_t(N,rho,recovery,1,ku,I,k,nu);
        p=@(y) tcdf((k-y*sqrt(rho*(nu-2)/nu))/sqrt((1-rho)*(nu-2)/nu),nu);
        density=@(y) tpdf(y,nu);
    case 1
        mezzanine=N*(1-ku)-KL_t(N,rho,recovery,1,ku,I,k,nu);
        p=@(eta) normcdf((eta./sqrt(1-rho)));
        integrand = @(eta) integral(@(w) exp(-1./(2*rho).*(eta-k.*sqrt(w/nu)).^2).*w.^(nu/2-1).*exp(-w/2),0,100, 'ArrayValued', true);
        density = @(eta) 1./(sqrt(rho*pi).*2^((nu+1)/2)*gamma(nu/2))*integrand(eta);
    case 2
        mezzanine=N*(1-ku)-KL_vasicek(N,rho,recovery,1,ku,I,k);
        p=@(y) normcdf((k-y*sqrt(rho))/sqrt(1-rho));
        density=@(y) normpdf(y);
end
lossrp=N*(1-recovery)*quadgk(@(y) p(y).*density(y), -100, +100);
price=N*ku-lossrp+mezzanine;

end