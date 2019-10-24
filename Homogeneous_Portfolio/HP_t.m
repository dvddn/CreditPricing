function price=HP_t(N,rho,recovery,ku,kd,I,k,nu)
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

lastwarn('')
warning('off', 'all')

u=ku/(1-recovery);
d=kd/(1-recovery);

% Default probability
p=@(eta) normcdf((eta/sqrt(1-rho))); 
% Loss function
L=@(z) min(max(z-d, 0), u-d)/(u-d);

exp_loss=0;
m=1;

while(m<=I && isempty(lastwarn))
    bin_coeff=nchoosek(I,m);
    integrand = @(eta) integral(@(w) exp(-1./(2*rho).*(eta-k.*sqrt(w/nu)).^2).*w.^(nu/2-1).*exp(-w/2),0,100,'ArrayValued', true);
    density = @(eta) 1./(sqrt(rho*pi).*2^((nu+1)/2)*gamma(nu/2))*integrand(eta);
    integrand2 = @(eta) bin_coeff*(p(eta).^m).*((1-p(eta)).^(I-m)).*density(eta);
    exp_loss= exp_loss+L(m/I).* integral(@(eta) arrayfun(integrand2, eta) ,-10,10, 'ArrayValued', true); 
    m=m+1;
end

warning('on', 'all')

if isempty(lastwarn)
    price=(ku-kd)*N*(1-exp_loss);
else
    price=[];
    warning('Too many obligors, price cannot be computed via exact formula')
end  

%In case the binomial gets too large, the function interrupts and shows a warning


end % function HP_t