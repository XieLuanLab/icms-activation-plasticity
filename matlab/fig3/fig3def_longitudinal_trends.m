%% Plot and count Fig 3D-F data points
% Reproduces the D-F scatter/trend plots and reports exact counts
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
% load(FIG3_DATA);  % uncomment if workspace not already loaded

metrics = {mergedDFF, mergedTrendsOnset, mergedTrendsDur};
metricNames = {'Maximum \DeltaF/F', 'Spike Onset (ms)', 'Spike Duration (ms)'};
panelNames = {'D', 'E', 'F'};

% Weekly timepoints (0-4)
if exist('weekTimepoints', 'var')
    wkTP = weekTimepoints;
else
    wkTP = 0:4;
end

for t = 1:3
    metricTrend = metrics{t};

    fprintf('\n========================================\n');
    fprintf('Panel %s: %s\n', panelNames{t}, metricNames{t});
    fprintf('========================================\n');
    fprintf('Total response units (rows in mergedDFF): %d\n', size(metricTrend, 1));

    %% Count valid data points per week
    fprintf('\nValid data points per week:\n');
    totalDataPoints = 0;
    for w = 1:5
        valid_w = isfinite(metricTrend(:,w)) & metricTrend(:,w) ~= 0;
        fprintf('  Week %d: %d data points\n', w-1, sum(valid_w));
        totalDataPoints = totalDataPoints + sum(valid_w);
    end
    fprintf('  Total data points plotted: %d\n', totalDataPoints);

    %% Count for early vs late comparison
    earlyData = mean(metricTrend(:, 1:2), 2, 'omitnan');
    lateData = mean(metricTrend(:, 3:5), 2, 'omitnan');
    valid = isfinite(earlyData) & isfinite(lateData);
    fprintf('\nEarly vs Late comparison:\n');
    fprintf('  Valid early+late: %d\n', sum(valid));
    fprintf('  Dropped (incomplete): %d\n', size(metricTrend,1) - sum(valid));

    %% Count unique lines (response units with any valid data)
    anyValid = any(isfinite(metricTrend) & metricTrend ~= 0, 2);
    fprintf('  Response units with any valid week: %d\n', sum(anyValid));

    %% Plot - reproduce D-F
    figure('Name', sprintf('Fig 3%s: %s', panelNames{t}, metricNames{t}), ...
        'Position', [100 + (t-1)*400, 200, 350, 300]);
    hold on

    % Plot each response unit
    metSkip = ~valid;
    for i = 1:size(metricTrend, 1)
        if ~metSkip(i)
            x = wkTP;
            y = metricTrend(i, :);
            validPts = isfinite(y) & y ~= 0;
            if sum(validPts) >= 2
                xv = x(validPts);
                yv = y(validPts);
                scatter1 = scatter(xv, yv, 'k');
                scatter1.MarkerFaceAlpha = 0.15;
                scatter1.MarkerEdgeAlpha = 0.15;
                P = polyfit(xv, yv, 1);
                yfit = P(1)*xv + P(2);
                plot(xv, yfit, 'color', [0 0 0 0.15], 'linewidth', 1);
            end
        end
    end

    % Overall trend line
    aveAll = mean(metricTrend(~metSkip, :), 1, 'omitnan');
    validAve = isfinite(aveAll);
    P = polyfit(wkTP(validAve), aveAll(validAve), 1);
    yfit = P(1)*wkTP + P(2);
    plot(wkTP, yfit, 'color', [0 0 0], 'linewidth', 3);

    xlim([min(wkTP)-0.5, max(wkTP)+0.5]);
    xticks(0:4);
    xlabel('Weeks of Training');
    ylabel(metricNames{t});
    title(sprintf('Fig 3%s (n=%d)', panelNames{t}, sum(valid)));
    box off

    %% Wilcoxon signed-rank test (paired)
    earlyValid = earlyData(valid);
    lateValid = lateData(valid);
    [p, ~, stats] = signrank(earlyValid, lateValid);
    r = abs(stats.zval) / sqrt(sum(valid));
    fprintf('\nWilcoxon signed-rank test:\n');
    fprintf('  n = %d, p = %.4g, r = %.4f\n', sum(valid), p, r);
end

%% Per-animal breakdown
fprintf('\n========================================\n');
fprintf('Per-animal breakdown of response units\n');
fprintf('========================================\n');
animalNames = {'ICMS92', 'ICMS98', 'ICMS93'};
for a = 1:3
    dff = allDFF{a};
    fprintf('%s: %d response units x %d sessions\n', animalNames{a}, size(dff, 1), size(dff, 2));

    % How many valid per weekly bin
    % Sessions -> weeks: need daysTrained to map
    days = allDays{a};
    weeks = ceil(days / 7);
    fprintf('  Session weeks: '); fprintf('%d ', weeks); fprintf('\n');
end

fprintf('\nTotal response units: %d\n', size(allDFF{1},1) + size(allDFF{2},1) + size(allDFF{3},1));
fprintf('After early+late filter: %d\n', sum(isfinite(mean(mergedDFF(:,1:2),2,'omitnan')) & isfinite(mean(mergedDFF(:,3:5),2,'omitnan'))));
