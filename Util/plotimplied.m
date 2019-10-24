function plotimplied(values, marketvalues, k)
% function plotting the comparative between market prices/implied corr and
% model prices/implied corr
% 
% INPUTS:
% values:       model implied correlations or prices
% marketvalues: market implied correlations or prices
% k:            detachment points

ax=gca;
plot(k,values, 'm^-', 'LineWidth', 1)
ax.XTick=[0.03,0.06,0.09,0.12,0.22];
ax.XLabel.String= 'Detachment points';
if max(values)<=1
    ax.YLabel.String= 'Implied correlation';
    ax.YLim=[0,0.5];
else
    ax.YLabel.String= 'Prices';
end
grid on
hold on
plot(k,marketvalues, 'bd-', 'LineWidth', 1)
legend('Model','Market');
hold off


end %function plotimplied