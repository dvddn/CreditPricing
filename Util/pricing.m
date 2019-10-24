function [exact,KL,LHP]=pricing(N,rho,recovery,ku,kd,p,nu,obligors,model)
% function pricing an MBS tranche for several number of obligors with three
% different methods: exact solution, KL approximation and LHP assumption.
% Pricing can be done according to three different models.
% 
% INPUTS:
% N:        MBS notional
% rho:      model correlation
% recovery: average recovery of each mortgage
% ku:       upper detachment point of the tranche
% kd:       lower detachment point of the tranche
% p:        default probability of each mortgage
% nu:       # of degrees of freedom of the model (unused if gaussian model)
% obligors: vector containing different numbers of obligors
% model:    =0 for double t, =1 for student t, =2 for vasicek
% 
% OUTPUTS:
% exact:    tranche prices computed with exact solution up to a computable number of obligors
% KL:       tranche prices computed with KL approx
% LHP:      tranche price computed with LHP assumption
% 
% USES:
% LHP_double_t
% HP_double_t
% KL_double_t
% LHP_t
% HP_t
% KL_t
% LHP_vasicek
% HP_vasicek
% KL_vasicek
% KL_equity

exact=[];
KL=zeros(length(obligors),1);

switch model
    case 0
        k=inversep(rho,p,nu);
        LHP=LHP_double_t(N,rho,recovery,ku,kd,k,nu);
        for j=1:length(obligors)
            exact= [exact;HP_double_t(N,rho,recovery,ku,kd,obligors(j),k,nu)];
            if kd
                KL(j)=KL_double_t(N,rho,recovery,ku,kd,obligors(j),k,nu);
            else
                KL(j)= KL_equity(N,rho,recovery,ku,obligors(j),k,model,nu);
            end
        end
        
    case 1
        k=tinv(p,nu);
        LHP=LHP_t(N,rho,recovery,ku,kd,k,nu);
        for j=1:length(obligors)
            exact= [exact;HP_t(N,rho,recovery,ku,kd,obligors(j),k,nu)];
            if kd
                KL(j)= KL_t(N,rho,recovery,ku,kd,obligors(j),k,nu);
            else
                KL(j)= KL_equity(N,rho,recovery,ku,obligors(j),k,model,nu);
            end
        end
        
    case 2
        k=norminv(p);
        LHP=LHP_vasicek(N,rho,recovery,ku,kd,k);
        for j=1:length(obligors)
            exact= [exact;HP_vasicek(N,rho,recovery,ku,kd,obligors(j),k)];
            if kd
                KL(j)= KL_vasicek(N,rho,recovery,ku,kd,obligors(j),k);
            else
                KL(j)= KL_equity(N,rho,recovery,ku,obligors(j),k,model);
            end
        end
end

end %function pricing