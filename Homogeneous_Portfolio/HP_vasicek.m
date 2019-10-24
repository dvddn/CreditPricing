function price=HP_vasicek(N,rho,recovery,ku,kd,I,k)
% Price of a tranche with the closed formula for Vasicek model under HP
% assumptions
%
%INPUT
% N:           Notional of the portfolio
% rho:         Correlation
% recovery:    Recovery rate
% ku:          Detachment point of the tranche
% kd:          Attachment point of the tranche
% I:           Number of mortgages
% k:           effective parameter (norminv(default probability of each mortgage))

lastwarn('');
warning('off', 'all')

u=ku/(1-recovery);
d=kd/(1-recovery);

% Default probability
p=@(y) normcdf((k-y*sqrt(rho))/sqrt(1-rho)); 
% Loss function
L=@(z) min(max(z-d, 0), u-d)/(u-d); 

exp_loss=0;
m=1;
% Evaluation of the expected loss via closed formula
while (m<=I && isempty(lastwarn))
    bin_coeff=nchoosek(I,m);
    integrand= @(y) normpdf(y).*bin_coeff.* (p(y).^m).* ((1-p(y)).^(I-m)); 
    exp_loss= exp_loss+L(m/I).* quadgk( integrand, -10, 10);
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

end % function HP_vasicek