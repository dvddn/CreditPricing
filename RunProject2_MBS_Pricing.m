%% MBS Pricing - André, Di Nello, Mini
% Project 2, AY 2016-17
clear all;
close all;
lastwarn('');
clc;

%% Initialization TBM
model = 0;  % 0: Double t    1: Single t
flag = 2;   % 1: Calibration through Compound Correlation
            % 2: Calibration through Base Correlation
            % 3: Calibration through Mezzanine Prices - just for double checking
            % 4: Calibration through Cumulated Equity Prices - just for double checking
tol=1e-4;   % calibration algorithm precision
maxIter=50;
int=[2+5*model,20-8*model]; % search interval [2 20] for double t
                             %                [7 12] for single t
N=2*10^9;   % notional amount
p=0.07;     % probability of default of each mortgage
recovery=0.65; % average recovery of every mortgage
[Data,Market] = initialization(N,recovery,p);
% inizialization of detachment points, market correlation and prices

%% Calibration
[nu] = Calibration(Data,Market,flag,model,tol,maxIter,int);
[rho, values, exitflag]= ImpliedComputation(Data,Market,model,flag,nu);

%% Pricing
obligors=[10 20 50 100 200 500 1000 2000 5000 10000];

for i=1:3
[exact, KL, LHP]=pricing(Data.N,rho,Data.recovery,Data.ku(i),Data.kd(i),Data.p,nu,obligors,model);
plot_results(obligors,exact,KL,LHP);
end

%% KL Calibration
Data.I=500
[nu2] = CalibrationKL(Data,Market,model,tol,maxIter,int);
[rho2, values2, exitflag2]= ImpliedComputationKL(Data,Market,model,nu2);
