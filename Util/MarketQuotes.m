function [Cumulative, Mezzanines, CompoundCorr] = MarketQuotes(N,Pi,ku,kd,k,MktCorr)
% Function computing Market Prices of Cumulative and Mezzanine tranches, and
% implied Compound Correlation
%
% INPUT:
% N:        Notional amount
% Pi:       average recovery of every mortgage
% ku:       detachment points of the mezzanine (and equity) tranches (upper notional bound)
% kd:       subordination points of the mezzanine tranches (lower notional bound)
% k:        effective parameter: notional of each credit
% 
% OUTPUTS:
% Cumulative:   Market value of the cumulated equity tranches
% Mezzanine:    Quoted Market value of mezzanine tranches
% CompundCorr:  Compound Implied correlation related to the mezzanine tranches
% 
% USES: 
% LHP_vasicek:  computes the price of a MBS tranche with gaussian model

if length(MktCorr)~=length(ku)&&length(MktCorr)~=length(kd)
    error('Market Correlations must correspond to detachment points')
end

Cumulative = zeros(length(MktCorr),1);
CompoundCorr=Cumulative;
CompoundCorr(1)=MktCorr(1);

Cumulative(1)=LHP_vasicek(N,MktCorr(1),Pi,ku(1),kd(1),k);
Mezzanines = Cumulative;

for i=2:length(MktCorr)
    Cumulative(i)=LHP_vasicek(N,MktCorr(i),Pi,ku(i),kd(1),k);
    Mezzanines(i)=Cumulative(i)-Cumulative(i-1);
    f=@(x) Mezzanines(i)-LHP_vasicek(N,x,Pi,ku(i),kd(i),k);
    try
        CompoundCorr(i)=fzero(f,MktCorr(i));
    catch
        CompoundCorr(i)=0;
    end
end


end % function MarketQuotes