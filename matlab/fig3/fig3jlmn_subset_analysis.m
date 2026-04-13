%% Plot and count Fig 3 J, L, M, N
% Requires fig3_final_results or fig3def_plot_data loaded
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end

if ~exist('mergedDFF', 'var')
    load(fullfile(DATA_ROOT, 'fig3def_plot_data.mat'));
end

subset_idx = subsetMerged == 1;
remaining_idx = subsetMerged == 0;

%% ========== Panel J: Subset vs remaining dF/F over weeks ==========
fprintf('\n========================================\n');
fprintf('Panel J: Subset vs remaining dF/F over weeks\n');
fprintf('========================================\n');

subDFF = mergedDFF(subset_idx, :);
remDFF = mergedDFF(remaining_idx, :);

fprintf('Subset: %d response units\n', size(subDFF, 1));
fprintf('Remaining: %d response units\n', size(remDFF, 1));

fprintf('\nPer-week valid counts (subset / remaining):\n');
for w = 1:5
    ns = sum(isfinite(subDFF(:,w)) & subDFF(:,w) ~= 0);
    nr = sum(isfinite(remDFF(:,w)) & remDFF(:,w) ~= 0);
    fprintf('  Week %d: %d / %d\n', w-1, ns, nr);
end

% Early vs late for each group (why are there responses with NaN entire
% early period)
earSub = mean(subDFF(:,1:2), 2, 'omitnan'); earSub = earSub(isfinite(earSub));
latSub = mean(subDFF(:,3:5), 2, 'omitnan'); latSub = latSub(isfinite(latSub));
earRem = mean(remDFF(:,1:2), 2, 'omitnan'); earRem = earRem(isfinite(earRem));
latRem = mean(remDFF(:,3:5), 2, 'omitnan'); latRem = latRem(isfinite(latRem));

fprintf('\nSubset early vs late:\n');
[p,~,stats] = ranksum(earSub, latSub);
r = abs(stats.zval)/sqrt(numel(earSub)+numel(latSub));
fprintf('  n_early=%d, n_late=%d, p=%.4g, r=%.4f\n', numel(earSub), numel(latSub), p, r);

fprintf('Remaining early vs late:\n');
[p,~,stats] = ranksum(earRem, latRem);
r = abs(stats.zval)/sqrt(numel(earRem)+numel(latRem));
fprintf('  n_early=%d, n_late=%d, p=%.4g, r=%.4f\n', numel(earRem), numel(latRem), p, r);

% Per-week subset vs remaining comparison
fprintf('\nPer-week subset vs remaining (Mann-Whitney U):\n');
for w = 1:5
    subWk = subDFF(:,w); subWk = subWk(isfinite(subWk) & subWk ~= 0);
    remWk = remDFF(:,w); remWk = remWk(isfinite(remWk) & remWk ~= 0);
    if ~isempty(subWk) && ~isempty(remWk)
        [p,~,stats] = ranksum(subWk, remWk);
        r = abs(stats.zval)/sqrt(numel(subWk)+numel(remWk));
        fprintf('  Week %d: n_sub=%d, n_rem=%d, p=%.4g, r=%.4f\n', w-1, numel(subWk), numel(remWk), p, r);
    end
end

% Plot
figure('Name', 'Fig 3J', 'Position', [100 200 500 350]);
hold on
col_sub = [0.8 0.2 0.2];
col_rem = [0.3 0.5 0.8];
wkTP = 0:4;

sub_mean = zeros(5,1); sub_std = zeros(5,1);
rem_mean = zeros(5,1); rem_std = zeros(5,1);
for w = 1:5
    sv = subDFF(:,w); sv = sv(isfinite(sv) & sv ~= 0);
    rv = remDFF(:,w); rv = rv(isfinite(rv) & rv ~= 0);
    if ~isempty(sv)
        sub_mean(w) = mean(sv); sub_std(w) = std(sv);
    end
    if ~isempty(rv)
        rem_mean(w) = mean(rv); rem_std(w) = std(rv);
    end
end

errorbar(wkTP-0.1, sub_mean, sub_std, '-o', 'Color', col_sub, 'LineWidth', 2, 'MarkerFaceColor', col_sub);
errorbar(wkTP+0.1, rem_mean, rem_std, '-o', 'Color', col_rem, 'LineWidth', 2, 'MarkerFaceColor', col_rem);
legend('Subset', 'Remaining', 'Box', 'off');
xlabel('Weeks of Training'); ylabel('\DeltaF/F_{max}');
title('Fig 3J'); xlim([-0.5 4.5]); xticks(0:4); box off;

%% ========== Panel L: Deconvolved spike count - subset vs remaining ==========
fprintf('\n========================================\n');
fprintf('Panel L: Deconvolved spike count\n');
fprintf('========================================\n');

subSpikes = mergedTrendsSpikes(subset_idx, :);
remSpikes = mergedTrendsSpikes(remaining_idx, :);

fprintf('Per-week valid counts (subset / remaining):\n');
for w = 1:5
    ns = sum(isfinite(subSpikes(:,w)) & subSpikes(:,w) ~= 0);
    nr = sum(isfinite(remSpikes(:,w)) & remSpikes(:,w) ~= 0);
    fprintf('  Week %d: %d / %d\n', w-1, ns, nr);
end

fprintf('\nPer-week subset vs remaining (Mann-Whitney U):\n');
for w = 1:5
    subWk = subSpikes(:,w); subWk = subWk(isfinite(subWk) & subWk ~= 0);
    remWk = remSpikes(:,w); remWk = remWk(isfinite(remWk) & remWk ~= 0);
    if ~isempty(subWk) && ~isempty(remWk)
        [p,~,stats] = ranksum(subWk, remWk);
        r = abs(stats.zval)/sqrt(numel(subWk)+numel(remWk));
        fprintf('  Week %d: n_sub=%d, n_rem=%d, p=%.4g, r=%.4f\n', w-1, numel(subWk), numel(remWk), p, r);
    end
end

% Plot
figure('Name', 'Fig 3L', 'Position', [550 200 500 350]);
hold on
sub_mean_sp = zeros(5,1); sub_std_sp = zeros(5,1);
rem_mean_sp = zeros(5,1); rem_std_sp = zeros(5,1);
for w = 1:5
    sv = subSpikes(:,w); sv = sv(isfinite(sv) & sv ~= 0);
    rv = remSpikes(:,w); rv = rv(isfinite(rv) & rv ~= 0);
    if ~isempty(sv)
        sub_mean_sp(w) = mean(sv); sub_std_sp(w) = std(sv);
    end
    if ~isempty(rv)
        rem_mean_sp(w) = mean(rv); rem_std_sp(w) = std(rv);
    end
end
errorbar(wkTP-0.1, sub_mean_sp, sub_std_sp, '-o', 'Color', col_sub, 'LineWidth', 2, 'MarkerFaceColor', col_sub);
errorbar(wkTP+0.1, rem_mean_sp, rem_std_sp, '-o', 'Color', col_rem, 'LineWidth', 2, 'MarkerFaceColor', col_rem);
legend('Subset', 'Remaining', 'Box', 'off');
xlabel('Weeks of Training'); ylabel('Deconvolved spikes');
title('Fig 3L'); xlim([-0.5 4.5]); xticks(0:4); box off;

%% ========== Panel M: Deconvolved spike onset - subset vs remaining ==========
fprintf('\n========================================\n');
fprintf('Panel M: Deconvolved spike onset\n');
fprintf('========================================\n');

subOnset = mergedTrendsSpikesItime(subset_idx, :);
remOnset = mergedTrendsSpikesItime(remaining_idx, :);

fprintf('Per-week valid counts (subset / remaining):\n');
for w = 1:5
    ns = sum(isfinite(subOnset(:,w)) & subOnset(:,w) ~= 0);
    nr = sum(isfinite(remOnset(:,w)) & remOnset(:,w) ~= 0);
    fprintf('  Week %d: %d / %d\n', w-1, ns, nr);
end

fprintf('\nPer-week subset vs remaining (Mann-Whitney U):\n');
for w = 1:5
    subWk = subOnset(:,w); subWk = subWk(isfinite(subWk) & subWk ~= 0);
    remWk = remOnset(:,w); remWk = remWk(isfinite(remWk) & remWk ~= 0);
    if ~isempty(subWk) && ~isempty(remWk)
        [p,~,stats] = ranksum(subWk, remWk);
        r = abs(stats.zval)/sqrt(numel(subWk)+numel(remWk));
        fprintf('  Week %d: n_sub=%d, n_rem=%d, p=%.4g, r=%.4f\n', w-1, numel(subWk), numel(remWk), p, r);
    end
end

% Plot
figure('Name', 'Fig 3M', 'Position', [1000 200 500 350]);
hold on
sub_mean_on = zeros(5,1); sub_std_on = zeros(5,1);
rem_mean_on = zeros(5,1); rem_std_on = zeros(5,1);
for w = 1:5
    sv = subOnset(:,w); sv = sv(isfinite(sv) & sv ~= 0);
    rv = remOnset(:,w); rv = rv(isfinite(rv) & rv ~= 0);
    if ~isempty(sv)
        sub_mean_on(w) = mean(sv); sub_std_on(w) = std(sv);
    end
    if ~isempty(rv)
        rem_mean_on(w) = mean(rv); rem_std_on(w) = std(rv);
    end
end
errorbar(wkTP-0.1, sub_mean_on, sub_std_on, '-o', 'Color', col_sub, 'LineWidth', 2, 'MarkerFaceColor', col_sub);
errorbar(wkTP+0.1, rem_mean_on, rem_std_on, '-o', 'Color', col_rem, 'LineWidth', 2, 'MarkerFaceColor', col_rem);
legend('Subset', 'Remaining', 'Box', 'off');
xlabel('Weeks of Training'); ylabel('Spike onset (ms)');
title('Fig 3M'); xlim([-0.5 4.5]); xticks(0:4); box off;

%% ========== Panel N: Bar plots - week0 subset vs weeks2-4 subset vs end-recruited ==========
fprintf('\n========================================\n');
fprintf('Panel N: Three-group comparison\n');
fprintf('========================================\n');

metricSets = {
    mergedDFF, mergedEndDFF, 'dF/Fmax';
    mergedTrendsSpikes, mergedEndTrendsSpikes, 'Spike count';
    mergedTrendsSpikesItime, mergedEndTrendsSpikesItime, 'Spike onset'
};

figure('Name', 'Fig 3N', 'Position', [100 100 900 300]);

for m = 1:size(metricSets, 1)
    metricAll = metricSets{m, 1};
    metricEnd = metricSets{m, 2};
    metricName = metricSets{m, 3};

    metricSub = metricAll(subset_idx, :);

    % Group 1: subset at week 0
    earSub = metricSub(:, 1);
    earSub = earSub(isfinite(earSub) & earSub ~= 0);

    % Group 2: subset at weeks 2-4 (averaged per neuron)
    latSub = mean(metricSub(:, 3:5), 2, 'omitnan');
    latSub = latSub(isfinite(latSub) & latSub ~= 0);

    % Group 3: end-recruited
    endVals = metricEnd(:);
    endVals = endVals(isfinite(endVals) & endVals ~= 0);

    fprintf('\n--- %s ---\n', metricName);
    fprintf('  Group 1 (subset week 0): n=%d, median=%.2f\n', numel(earSub), median(earSub));
    fprintf('  Group 2 (subset weeks 2-4): n=%d, median=%.2f\n', numel(latSub), median(latSub));
    fprintf('  Group 3 (end-recruited): n=%d, median=%.2f\n', numel(endVals), median(endVals));

    % Kruskal-Wallis
    barDat = [earSub; latSub; endVals];
    barG = [ones(numel(earSub),1); 2*ones(numel(latSub),1); 3*ones(numel(endVals),1)];
    [p_kw, ~, stats_kw] = kruskalwallis(barDat, barG, 'off');
    c = multcompare(stats_kw, 'Display', 'off');

    fprintf('  Kruskal-Wallis: p=%.4g\n', p_kw);
    grpNames = {'wk0_subset', 'wk2-4_subset', 'end_recruited'};
    for k = 1:size(c, 1)
        fprintf('  %s vs %s: p=%.4g\n', grpNames{c(k,1)}, grpNames{c(k,2)}, c(k,6));
    end

    % Bar plot
    subplot(1, 3, m);
    data = [mean(earSub), mean(latSub), mean(endVals)];
    stds = [std(earSub), std(latSub), std(endVals)];
    bar(1:3, data);
    hold on
    errorbar(1:3, data, stds, 'k', 'LineStyle', 'none');
    set(gca, 'XTickLabel', {'Wk0 sub', 'Wk2-4 sub', 'End rec'});
    ylabel(metricName);
    title(sprintf('%s (KW p=%.2g)', metricName, p_kw));
    box off;
end

%% ========== FULL SUMMARY ==========
fprintf('\n========================================\n');
fprintf('FULL SUMMARY OF COUNTS\n');
fprintf('========================================\n');
fprintf('Panel I:\n');
fprintf('  First valid dF/F histogram: %d data points\n', sum(valid_bMax));
fprintf('  Last valid dF/F histogram: %d data points (%d subset, %d remaining)\n', ...
    sum(valid_fMax), sum(valid_fMax & subset_idx), sum(valid_fMax & remaining_idx));

fprintf('Panel J:\n');
fprintf('  Subset: %d total (%d early, %d late)\n', sum(subset_idx), numel(earSub), numel(latSub));
fprintf('  Remaining: %d total (%d early, %d late)\n', sum(remaining_idx), numel(earRem), numel(latRem));

fprintf('Panel L-M: same groups as J, using spike count and spike onset\n');

fprintf('Panel N:\n');
fprintf('  Group 1 (subset week 0): %d\n', numel(earSub));
fprintf('  Group 2 (subset weeks 2-4): %d\n', numel(latSub));
fprintf('  Group 3 (end-recruited): %d\n', numel(endVals));
