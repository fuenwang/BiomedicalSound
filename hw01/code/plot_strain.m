fig = figure();
subplot('311')
plot(depth_10_0(1:end-1)*1e3, strain_10_0 * 100, 'Linewidth', 2)
title('Thermal Strain(%)(M=10, N=0)')
xlabel('Depth(mm)', 'FontSize', 15)
ylabel('\boldmath{$\frac{\partial \Delta t(z)}{\partial z}$ (\%)}', 'interpreter', 'latex', 'FontSize', 15)

subplot('312')
plot(depth_10_50(1:end-1)*1e3, strain_10_50 * 100, 'Linewidth', 2)
title('Thermal Strain(%)(M=10, N=0.5)')
xlabel('Depth(mm)', 'FontSize', 15)
ylabel('\boldmath{$\frac{\partial \Delta t(z)}{\partial z}$ (\%)}', 'interpreter', 'latex', 'FontSize', 15)

subplot('313')
plot(depth_10_75(1:end-1)*1e3, strain_10_75 * 100, 'Linewidth', 2)
title('Thermal Strain(%)(M=10, N=0.75)')
xlabel('Depth(mm)', 'FontSize', 15)
ylabel('\boldmath{$\frac{\partial \Delta t(z)}{\partial z}$ (\%)}', 'interpreter', 'latex', 'FontSize', 15)

fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../doc/src/strain_10.pdf', '-dpdf')