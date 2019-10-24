function plot_results(obligors,exact,KL,LHP)
% function plotting results of function "pricing"
% INPUTS:
% obligors: vector containing different numbers of obligors
% exact:    vector of exact prices for each number of obligors (if possible)
% KL:       vector of KL approx prices for each number of obligors
% LHP:      scalar LHP price

figure
ax=gca;
semilogx(obligors(1:length(exact)),exact, 'b-d', 'LineWidth', 1.0)
hold on
grid on
semilogx(obligors, KL, 'r-^', 'LineWidth',1.0)
semilogx([10 1e4], [1 1]*LHP, 'g', 'LineWidth',1.0)
legend('Exact solution', 'KL approx solution', 'LHP approx solution','Location','East')
ax.XLabel.String= 'Number of obligors';
ax.YLabel.String= 'Price';
hold off

end %function plot_results