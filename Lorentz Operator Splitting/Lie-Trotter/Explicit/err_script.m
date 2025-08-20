
addpath('../..')

K = 1;

[delx_full,dely_full,delt_full,X_full,Y_full,H_full] = LT_2d_Lorentz(32*K);

%%

N = [1,2,4,8].*K;

delx = zeros(size(N));
dely = zeros(size(N));
delt = zeros(size(N));
X = {};
Y = {};
H = {};

guide_line_1 = 1./N;
guide_line_2 = 1./(N.^2);

err = zeros(size(N));

for i = 1:length(N)

    [delx(i),dely(i),delt(i),X{i},Y{i},H{i}] = LT_2d_Lorentz(N(i));
end
%%

for i = 1:length(N)
    interp_H = interp2(X_full,Y_full,H_full,X{i},Y{i},'spline');
    err(i) = (norm(interp_H - H{i}))*sqrt(X{i}(end,2) - X{i}(end,1))...
        *sqrt(Y{i}(2,end) - Y{i}(1,end));
end

%%

% figure(3)
% hold on
% plot(x_full,E_full)
% title('Final Snap Shot Plot')
% xlabel('X')
% ylabel('Electric Field')
% legend('1','2','3','4','full')
% hold off


%%
line_diff_1 = log(err(1)) - log(guide_line_1(1))+ 0.5;
line_diff_2 = log(err(1)) - log(guide_line_2(1))- 0.5;

figure(6)
plot(log(delx),log(err),'-o')
hold on
plot(log(delx),log(guide_line_1) + line_diff_1,':',Color="k")
plot(log(delx),log(guide_line_2) + line_diff_2,'--',Color="k")
title('Log-Log Plot of Final Snapshot Error')
xlabel('log(\Delta x)')
ylabel('log(L_2 Error)')
legend('Numerical Error','O(h) guideline','O(h^2) guideline','Location','southeast')
hold off
disp("Final Snapshot order of convergence")
disp((log(err(end)/err(end-1)))/(log(delx(end)/delx(end-1))))
savefig('err_plot_grid_refine.fig')
saveas(gcf,'err_plot_grid_refine.png')

%%

out_matrix = zeros(length(N),6);
out_matrix(:,1) = delx;
out_matrix(:,2) = dely;
out_matrix(:,3) = delt;

out_matrix(:,4) = err;
out_matrix(2:end,5) = err(1:end-1)./err(2:end);
out_matrix(2:end,6) = log2(out_matrix(2:end,5));
out_matrix(1,5) = nan;
out_matrix(1,6) = nan;

matrix2latex(out_matrix,'err_table_grid_refine.tex','columnLabels', {'$\Delta x$', '$\Delta y$', '$\Delta t$', '$L_2$ Error', 'Ratio', 'Rate'})
