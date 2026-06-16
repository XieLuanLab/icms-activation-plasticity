%% Fig3 Statistics Extraction
% Loads Plotter2Results and computes all stats for Fig3 panels
% Outputs to fig3/fig3_stats.csv
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
load(FIG3_DATA);

csvPath = fullfile(fileparts(mfilename('fullpath')), 'fig3_stats.csv');
if exist(csvPath, 'file'), delete(csvPath); end

% Write header
fid = fopen(csvPath, 'w');
fprintf(fid, 'panel,comparison,test,n1,n2,p,effect_size_r\n');
fclose(fid);

%% Fig 3D-F: DFF, Onset, Duration - early vs late (Mann-Whitney U)
metrics = {mergedDFF, mergedTrendsOnset, mergedTrendsDur};
names = {'Fig3D', 'Fig3E', 'Fig3F'};
labels = {'dF/Fmax', 'onset_time', 'spike_duration'};

for t = 1:3
    metricTrend = metrics{t};
    % Per-neuron mean of early (weeks 0-1) vs late (weeks 2-4), matching roiPLotter2.m L4111
    earlyData = mean(metricTrend(:, 1:2), 2, 'omitnan');
    lateData = mean(metricTrend(:, 3:5), 2, 'omitnan');
    valid = isfinite(earlyData) & isfinite(lateData);
    earlyData = earlyData(valid);
    lateData = lateData(valid);

    [p, ~, stats] = ranksum(earlyData, lateData);
    n1 = numel(earlyData); n2 = numel(lateData);
    if isfield(stats, 'zval')
        r = abs(stats.zval) / sqrt(n1 + n2);
    else
        r = NaN;
    end

    fid = fopen(csvPath, 'a');
    fprintf(fid, '%s,%s early(0-1) vs late(2-4),Mann-Whitney U,%d,%d,%.4g,%.4f\n', ...
        names{t}, labels{t}, n1, n2, p, r);
    fclose(fid);
end

%% Fig 3J: Subset vs remaining - DFF early vs late per group
subset_idx = subsetMerged == 1;
remaining_idx = subsetMerged == 0;

% Subset (high-intensity) early vs late
subDFF = mergedDFF(subset_idx, :);
earSub = mean(subDFF(:, 1:2), 2, 'omitnan'); earSub = earSub(isfinite(earSub));
latSub = mean(subDFF(:, 3:5), 2, 'omitnan'); latSub = latSub(isfinite(latSub));
[p, ~, stats] = ranksum(earSub, latSub);
n1 = numel(earSub); n2 = numel(latSub);
if isfield(stats, 'zval'), r = abs(stats.zval)/sqrt(n1+n2); else, r = NaN; end
fid = fopen(csvPath, 'a');
fprintf(fid, 'Fig3J,subset_dF/F early vs late,Mann-Whitney U,%d,%d,%.4g,%.4f\n', n1, n2, p, r);
fclose(fid);

% Remaining early vs late
remDFF = mergedDFF(remaining_idx, :);
earRem = mean(remDFF(:, 1:2), 2, 'omitnan'); earRem = earRem(isfinite(earRem));
latRem = mean(remDFF(:, 3:5), 2, 'omitnan'); latRem = latRem(isfinite(latRem));
[p, ~, stats] = ranksum(earRem, latRem);
n1 = numel(earRem); n2 = numel(latRem);
if isfield(stats, 'zval'), r = abs(stats.zval)/sqrt(n1+n2); else, r = NaN; end
fid = fopen(csvPath, 'a');
fprintf(fid, 'Fig3J,remaining_dF/F early vs late,Mann-Whitney U,%d,%d,%.4g,%.4f\n', n1, n2, p, r);
fclose(fid);

% Per-week subset vs remaining
for wk = 1:5
    subWk = subDFF(:, wk); subWk = subWk(isfinite(subWk) & subWk ~= 0);
    remWk = remDFF(:, wk); remWk = remWk(isfinite(remWk) & remWk ~= 0);
    if ~isempty(subWk) && ~isempty(remWk)
        [p, ~, stats] = ranksum(subWk, remWk);
        n1 = numel(subWk); n2 = numel(remWk);
        if isfield(stats, 'zval'), r = abs(stats.zval)/sqrt(n1+n2); else, r = NaN; end
        fid = fopen(csvPath, 'a');
        fprintf(fid, 'Fig3J,subset vs remaining week%d dF/F,Mann-Whitney U,%d,%d,%.4g,%.4f\n', wk-1, n1, n2, p, r);
        fclose(fid);
    end
end

%% Fig 3L: Deconvolved spike count - per-week subset vs remaining
subSpikes = mergedTrendsSpikes(subset_idx, :);
remSpikes = mergedTrendsSpikes(remaining_idx, :);
for wk = 1:5
    subWk = subSpikes(:, wk); subWk = subWk(isfinite(subWk) & subWk ~= 0);
    remWk = remSpikes(:, wk); remWk = remWk(isfinite(remWk) & remWk ~= 0);
    if ~isempty(subWk) && ~isempty(remWk)
        [p, ~, stats] = ranksum(subWk, remWk);
        n1 = numel(subWk); n2 = numel(remWk);
        if isfield(stats, 'zval'), r = abs(stats.zval)/sqrt(n1+n2); else, r = NaN; end
        fid = fopen(csvPath, 'a');
        fprintf(fid, 'Fig3L,subset vs remaining week%d spike_count,Mann-Whitney U,%d,%d,%.4g,%.4f\n', wk-1, n1, n2, p, r);
        fclose(fid);
    end
end

%% Fig 3M: Deconvolved spike onset - per-week subset vs remaining
subOnset = mergedTrendsSpikesItime(subset_idx, :);
remOnset = mergedTrendsSpikesItime(remaining_idx, :);
for wk = 1:5
    subWk = subOnset(:, wk); subWk = subWk(isfinite(subWk) & subWk ~= 0);
    remWk = remOnset(:, wk); remWk = remWk(isfinite(remWk) & remWk ~= 0);
    if ~isempty(subWk) && ~isempty(remWk)
        [p, ~, stats] = ranksum(subWk, remWk);
        n1 = numel(subWk); n2 = numel(remWk);
        if isfield(stats, 'zval'), r = abs(stats.zval)/sqrt(n1+n2); else, r = NaN; end
        fid = fopen(csvPath, 'a');
        fprintf(fid, 'Fig3M,subset vs remaining week%d spike_onset,Mann-Whitney U,%d,%d,%.4g,%.4f\n', wk-1, n1, n2, p, r);
        fclose(fid);
    end
end

%% Fig 3N: Kruskal-Wallis comparisons (3 groups: subset week0, subset weeks2-4, end-recruited)
% Matches roiPLotter2.m group definitions (lines 4596-4603)
metricSets = {
    mergedDFF, mergedEndDFF, 'dF/Fmax';
    mergedTrendsSpikes, mergedEndTrendsSpikes, 'spike_count';
    mergedTrendsSpikesItime, mergedEndTrendsSpikesItime, 'spike_onset'
};

for m = 1:size(metricSets, 1)
    metricAll = metricSets{m, 1};
    metricEnd = metricSets{m, 2};
    metricName = metricSets{m, 3};

    metricSub = metricAll(subset_idx, :);

    % Group 1: SUBSET neurons at week 0 (matching roiPLotter2.m)
    earSub = metricSub(:, 1); earSub = earSub(isfinite(earSub) & earSub ~= 0);
    % Group 2: subset at weeks 2-4 (averaged per neuron)
    latSub = mean(metricSub(:, 3:5), 2, 'omitnan'); latSub = latSub(isfinite(latSub) & latSub ~= 0);
    % Group 3: end-recruited neurons
    endVals = metricEnd(:); endVals = endVals(isfinite(endVals) & endVals ~= 0);

    barDat = [earSub; latSub; endVals];
    barG = [ones(numel(earSub),1); 2*ones(numel(latSub),1); 3*ones(numel(endVals),1)];

    [p_kw, ~, stats_kw] = kruskalwallis(barDat, barG, 'off');
    c = multcompare(stats_kw, 'Display', 'off');

    fid = fopen(csvPath, 'a');
    fprintf(fid, 'Fig3N,%s overall Kruskal-Wallis,Kruskal-Wallis,%d,%d,%.4g,\n', ...
        metricName, numel(barDat), 3, p_kw);

    % Post-hoc pairwise from multcompare: c columns = [g1 g2 lower diff upper p]
    grpNames = {'week0_subset', 'weeks2-4_subset', 'end_recruited'};
    for k = 1:size(c, 1)
        g1 = c(k, 1); g2 = c(k, 2); p_post = c(k, 6);
        n1_post = sum(barG == g1); n2_post = sum(barG == g2);
        fprintf(fid, 'Fig3N,%s %s vs %s,Tukey-Kramer post-hoc,%d,%d,%.4g,\n', ...
            metricName, grpNames{g1}, grpNames{g2}, n1_post, n2_post, p_post);
    end
    fclose(fid);
end

fprintf('\nFig3 stats saved to: %s\n', csvPath);
