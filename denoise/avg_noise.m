function [mre_median_alpha, mre_alpha, y] = avg_noise
n = 128;
num = n * n;
error = importdata('denoise_error_pr01.mat');


fprintf('============ Euclidean\n');
mre = sum(error.error_euclid)/sum(error.num_euclid);
y(1) = mre;


fprintf('============ Log-Euclidean\n');
mre = sum(error.error_logeuclid)/sum(error.num_logeuclid);
y(2) = mre;


fprintf('============ symmldbreg\n');
mre = sum(error.error_symmldbreg)/sum(error.num_symmldbreg);
y(3) = mre;


fprintf('============ riemannian\n');
mre = sum(error.error_riemann)/sum(error.num_riemann);
y(4) = mre;


fprintf('============ riemannianmedian\n');
mre_riemannianmedian = sum(error.error_medianriemann)/sum(error.num_medianriemann);
y(5) = mre_riemannianmedian;


mre_alpha(1) = sum(error.error_alpha0)/sum(error.num_alpha0);
mre_alpha(2) = sum(error.error_alpha1)/sum(error.num_alpha1);
mre_alpha(3) = sum(error.error_alpha2)/sum(error.num_alpha2);
mre_alpha(4) = sum(error.error_alpha3)/sum(error.num_alpha3);
mre_alpha(5) = sum(error.error_alpha4)/sum(error.num_alpha4);
mre_alpha(6) = sum(error.error_alpha5)/sum(error.num_alpha5);
mre_alpha(7) = sum(error.error_alpha6)/sum(error.num_alpha6);
mre_alpha(8) = sum(error.error_alpha7)/sum(error.num_alpha7);
mre_alpha(9) = sum(error.error_alpha8)/sum(error.num_alpha8);
mre_alpha(10) = sum(error.error_alpha9)/sum(error.num_alpha9);

y(6) = mre_alpha(1);


mre_median_alpha(1) = sum(error.error_medianalpha0)/sum(error.num_medianalpha0);
mre_median_alpha(2) = sum(error.error_medianalpha1)/sum(error.num_medianalpha1);
mre_median_alpha(3) = sum(error.error_medianalpha2)/sum(error.num_medianalpha2);
mre_median_alpha(4) = sum(error.error_medianalpha3)/sum(error.num_medianalpha3);
mre_median_alpha(5) = sum(error.error_medianalpha4)/sum(error.num_medianalpha4);
mre_median_alpha(6) = sum(error.error_medianalpha5)/sum(error.num_medianalpha5);
mre_median_alpha(7) = sum(error.error_medianalpha6)/sum(error.num_medianalpha6);
mre_median_alpha(8) = sum(error.error_medianalpha7)/sum(error.num_medianalpha7);
mre_median_alpha(9) = sum(error.error_medianalpha8)/sum(error.num_medianalpha8);
mre_median_alpha(10) = sum(error.error_medianalpha9)/sum(error.num_medianalpha9);

y(7) = mre_median_alpha(2);





figure(1);
color_bar = lines(7);
y_total = y;
h = bar(1,y_total(1));
h.FaceColor = [0.6, 0.6, 0.6];
hold on;
for i = 2 : 7
    h = bar(i,y_total(i));
    h.FaceColor = color_bar(i-1, :);
end

%axis([0 11 0 45])
for i = 1 : 7
text(i,y_total(i), [num2str(y_total(i),'%0.3f')],...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom','FontSize',20,'Interpreter','Latex')
end
fname = {'Euc.','Log-Euc.', 'J-Div.','Rie.', 'Rie.-Med.', '$\alpha$-Div.', '$\alpha$-Med.'}; 
set(gca,'TickLength',[0 0], ...
    'XTick', 1:7, 'XTickLabel',fname, ...
    'FontSize',25,'ygrid','on', 'TickLabelInterpreter', 'Latex');
xlabel('Means','FontSize',25,'Interpreter','Latex');
ylabel('MRE','FontSize',25,'Interpreter','Latex');
title('Pr = 0.1','FontSize',25,'Interpreter','Latex')
hold off;
end