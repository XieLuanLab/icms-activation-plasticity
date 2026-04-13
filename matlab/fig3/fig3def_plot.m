%% Fig 3D-F - Plot longitudinal trends from saved data
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
load(fullfile(DATA_ROOT, 'fig3def_plot_data.mat'));

% Apply current filter (already NaN in data, just use as-is)
metrics = {mergedDFF, mergedTrendsOnset, mergedTrendsDur};
metricNames = {'Maximum dF/F', 'Spike Onset (ms)', 'Spike Duration (ms)'};
panelNames = {'D', 'E', 'F'};
wkTP = 0:4;

for t = 1:3
    m = metrics{t};
    early = mean(m(:,1:2), 2, 'omitnan');
    late = mean(m(:,3:5), 2, 'omitnan');
    valid = isfinite(early) & isfinite(late);

    figure('Name', sprintf('Fig 3%s', panelNames{t}), 'Position', [100+(t-1)*400 200 350 300]);
    hold on
    for i = 1:size(m,1)
        if ~valid(i), continue; end
        x = wkTP; y = m(i,:);
        vp = isfinite(y) & y ~= 0;
        if sum(vp) >= 2
            s1 = scatter(x(vp), y(vp), 'k');
            s1.MarkerFaceAlpha = 0.15; s1.MarkerEdgeAlpha = 0.15;
            P = polyfit(x(vp), y(vp), 1);
            plot(x(vp), P(1)*x(vp)+P(2), 'color', [0 0 0 0.15], 'linewidth', 1);
        end
    end
    aveAll = mean(m(valid,:), 1, 'omitnan');
    va = isfinite(aveAll);
    P = polyfit(wkTP(va), aveAll(va), 1);
    plot(wkTP, P(1)*wkTP+P(2), 'color', [0 0 0], 'linewidth', 3);
    xlim([-0.5 4.5]); xticks(0:4);
    xlabel('Weeks of Training'); ylabel(metricNames{t});
    title(sprintf('Fig 3%s (n=%d)', panelNames{t}, sum(valid)));
    box off

    % Wilcoxon signed-rank
    [p,~,stats] = signrank(early(valid), late(valid));
    r = abs(stats.zval)/sqrt(sum(valid));
    fprintf('%s: n=%d, p=%.4g, r=%.2f\n', panelNames{t}, sum(valid), p, r);
end
