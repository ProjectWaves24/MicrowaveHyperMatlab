% run first the program Precomputed_Inverse_Scattering_SOlution_M.m
close all

set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');

% x-axis to plot singular varlues

x = linspace(0, 255, 256);

% plotting singular values
plot(x', diag(SS),'LineWidth',2);

% L2-norm of right hand side d
L2d = 0.01;

legend('singular values  of \Sigma')


figure

% estimate for lambda=0

  bound= (L2d/diag(SS));

plot(x', bound,'LineWidth',2);

legend('estimate for \lambda=0')

font_size = 10;
set(gca, "FontSize", font_size)

set(gcf, "Units", "Inches", "Position", [0, 0, 7, 7], ...
       "PaperUnits", "Inches", "PaperSize", [7, 7])
