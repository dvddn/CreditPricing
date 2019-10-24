function [k]=inversep(rho,p,nu)
% function evaluating the default threshold k for the double t model
% INPUT:
% rho: model correlation
% p: default probability of each obligor
% nu: degrees of freedom of the double t model
% 
% OUTPUT:
% k: default threshold of the double t model

k2=tinv(p,nu); %used as a reference value for the search
f= @(y) quadgk(@(z)tcdf((y-sqrt(rho*(nu-2)/nu)*z)/sqrt((1-rho)*(nu-2)/nu),nu).*tpdf(z,nu),-100,100)-p;
k=fzero(f,k2);

end %function inversep