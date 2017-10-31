fig = figure();
subplot('311')
plot(depth_6_0(1:end-1)*1e3, strain_6_0 * 20, 'Linewidth', 2)
title('Thermal Strain(%)(M=6, N=0)')
xlabel('Depth(mm)', 'FontSize', 15)
ylabel('\boldmath{$\frac{\partial \Delta t(z)}{\partial z}$ (\%)}', 'interpreter', 'latex', 'FontSize', 15)

subplot('312')
plot(depth_6_50(1:end-1)*1e3, strain_6_50 * 20, 'Linewidth', 2)
title('Thermal Strain(%)(M=6, N=0.5)')
xlabel('Depth(mm)', 'FontSize', 15)
ylabel('\boldmath{$\frac{\partial \Delta t(z)}{\partial z}$ (\%)}', 'interpreter', 'latex', 'FontSize', 15)

subplot('313')
plot(depth_6_75(1:end-1)*1e3, strain_6_75 * 20, 'Linewidth', 2)
title('Thermal Strain(%)(M=6, N=0.75)')
xlabel('Depth(mm)', 'FontSize', 15)
ylabel('\boldmath{$\frac{\partial \Delta t(z)}{\partial z}$ (\%)}', 'interpreter', 'latex', 'FontSize', 15)

fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../doc/src/strain_6.pdf', '-dpdf')