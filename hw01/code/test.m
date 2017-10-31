fig = figure();
subplot('311')
plot(depth_10_0 * 1e3, delay_10_0 / 2 * 1e6, 'Linewidth', 2)
title('Echo-Time Shift(M=10, N=0%)')
xlabel('Depth(mm)')
ylabel('Delay(us)')

subplot('312')
plot(depth_10_50 * 1e3, delay_10_50 / 2 * 1e6, 'Linewidth', 2)
title('Echo-Time Shift(M=10, N=50%)')
xlabel('Depth(mm)')
ylabel('Delay(us)')

subplot('313')
plot(depth_10_75 * 1e3, delay_10_75 / 2 * 1e6, 'Linewidth', 2)
title('Echo-Time Shift(M=10, N=75%)')
xlabel('Depth(mm)')
ylabel('Delay(us)')

fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, '../doc/src/shift_10.pdf', '-dpdf')