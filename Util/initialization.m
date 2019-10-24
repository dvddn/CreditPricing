function [Data,Market]=initialization(N,recovery,p)
% function that initializes the fixed imputs of the problem
%
% INPUTS:
% N: Notional of the contract
% recovery: average recovery of each mortgage
% p: probability of default of each mortgage
% 
% OUTPUTS:
% Data: Struct containing input data +
%   ku: upper detachment points of the MBS tranches
%   kd: lower detachment points of the MBS tranches
%   k1: default threshold for the Vasicek model
%   I: default number of obligors
% Market: Struct containing:
%   base: given base market implied correlation
%   comp: calculated compound market implied correlation
%   eq: calculated market prices of the cumulated equity tranches
%   mezz: quoted market prices of the mezzanine tranches
% 
% USES:
% MarketQuotes: calculates Equity and Mezzanine tranche values, Compound
%               Correlation, starting from Market base correlation
% MktCorr.mat: Mat file containing market base correlations

ku=[0.03; 0.06; 0.09; 0.12; 0.22]; % standard values for upper detachment points
kd=[0; 0.03; 0.06; 0.09; 0.12];% standard values for lower detachment points
k1=norminv(p); 
temp=load('MktCorr.mat');
MktCorr=temp.MktCorr;
[Equity, Mezzanines, CompoundCorr] = MarketQuotes(N,recovery,ku,kd,k1,MktCorr);

Market=struct;
Market.base = MktCorr;
Market.comp = CompoundCorr;
Market.eq = Equity;
Market.mezz = Mezzanines;

Data=struct;
Data.N=N;
Data.p=p;
Data.recovery=recovery;
Data.k1=k1;
Data.ku=ku;
Data.kd=kd;
Data.I=1000;

end %function initialization