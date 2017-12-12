%
% EE6265 Fu-En Wang 106061531 HW2 11/14/2017
%

function saveFig(fig, path)
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig, path, '-dpdf')
end