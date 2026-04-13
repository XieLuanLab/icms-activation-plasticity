%% Fig 2G (threshold-only)
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
load(FIG2_DATA);

fontsize = 8;
set(groot,'DefaultTextInterpreter','tex', ...
          'DefaultAxesTickLabelInterpreter','tex', ...
          'DefaultLegendInterpreter','tex', ...
          'DefaultTextFontName','Arial', ...
          'DefaultAxesFontName','Arial', ...
          'DefaultLegendFontName','Arial', ...
          'DefaultAxesFontSize', fontsize, ...
          'DefaultTextFontSize', fontsize, ...
          'DefaultLegendFontSize', fontsize, ...
          'DefaultColorbarFontSize', fontsize);

blue = [0.298, 0.447, 0.690];

restrictThresholds = false;     % accept any finite threshold
validRange         = [4 6];    % µA, inclusive (unused when restrictThresholds = false)
dDiv               = 7;        % days per week bin

% Build the week grid (0..4) from all animals' days
daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
weekGrid           = unique(ceil(daysTrainedOverall / dDiv));
weekGrid           = weekGrid(weekGrid >= 0 & weekGrid <= 4);   

maxTraces = numel(anID) * max(cellfun(@numel, uChanAll));
mergingTraces = NaN(numel(weekGrid), maxTraces);
traceNum = 0;

figure(); hold on
ax = gca;

for aInd = 1:numel(anID)
    uChan          = uChanAll{aInd};          % [channels]
    uCurr          = uCurrAll{aInd};          % [currents available], e.g., [2 4 5 6 7 10]
    groupedCounts  = CurrChansAll{aInd};      % [days x channels x currentsIndex]
    daysTrained    = daysTrainedAll{aInd}(:); % [days]
    daysThresholds = daysThresholdsAll{aInd}; % [channels x days], µA

    uChanCurrs = uChanCurrsRef{aInd};
    popDates = popDatesAll{aInd};
    populationDistances = populationDistancesAll{aInd};

    if isempty(uChan) || isempty(uCurr) || isempty(groupedCounts), continue; end

    for chInd = 1:numel(uChan)
        % collect (week, value) for this channel using day-specific thresholds
        wk_vals = []; wk_bins = [];

        for d = 1:numel(daysTrained)
            curThresh = round(daysThresholds(chInd, d));

            % validity checks
            if ~isfinite(curThresh), continue; end
            if restrictThresholds && (curThresh < validRange(1) || curThresh > validRange(2))
                continue;
            end

            % find current index that matches today's threshold
            curIdx = find(uCurr == curThresh, 1, 'first');
            if isempty(curIdx), continue; end

            val = NaN;
            pIdx = find(uChanCurrs(:,1) == uChan(chInd) & uChanCurrs(:,2) == curThresh, 1, 'first');
            if ~isempty(pIdx)
                nDates = popDates{pIdx};
                dInd = find(nDates == d, 1, 'first');
                if ~isempty(dInd)
                    distances = populationDistances{pIdx}{dInd};
                    if ~isempty(distances)
                        distances = distances(:);
                        distances = distances(isfinite(distances) & distances > 0 & distances <= 800);
                        val = numel(distances);
                    end
                end
            end

            if ~isfinite(val) || val <= 0, continue; end

            wk = ceil(daysTrained(d) / dDiv);
            if wk < 0 || wk > 4, continue; end

            wk_vals(end+1,1) = val; 
            wk_bins(end+1,1) = wk; 
        end

        if isempty(wk_vals), continue; end

        % per-week median for this channel
        uW = unique(wk_bins,'stable');
        ch_med = zeros(numel(uW),1);
        for i = 1:numel(uW)
            ch_med(i) = median(wk_vals(wk_bins == uW(i)), 'omitnan');
        end

        % scatter the channel's per-week medians
        scatter(uW, ch_med, 15, blue, 'filled', 'MarkerFaceAlpha', 0.7);

        % align into merge matrix
        alignedRows = nan(numel(uW),1);
        for i = 1:numel(uW)
            r = find(weekGrid == uW(i), 1, 'first');
            if ~isempty(r), alignedRows(i) = r; end
        end
        use = isfinite(alignedRows);
        if any(use)
            traceNum = traceNum + 1;
            mergingTraces(alignedRows(use), traceNum) = ch_med(use);
        end
    end
end

% Trim unused columns
mergingTraces = mergingTraces(:, 1:traceNum);

% Per-week median across traces
curAve = median(mergingTraces, 2, 'omitnan');

% Plot median line
plot(weekGrid, curAve, '-', 'Color', blue, 'LineWidth', 2);

xlabel('Weeks of training');
ylabel('Number of neurons');
title('At threshold', 'FontWeight', 'normal');
xlim([-0.5, 4.5]); xticks(0:4);
ylim([0, max([curAve(:); mergingTraces(:)], [], 'omitnan') * 1.15]);
box off

% Early (0–1) vs Late (2–4) stats
early_rows = (weekGrid >= 0 & weekGrid <= 1);
late_rows  = (weekGrid >= 2 & weekGrid <= 4);
early_vals = mergingTraces(early_rows, :); early_vals = early_vals(:); early_vals = early_vals(isfinite(early_vals));
late_vals  = mergingTraces(late_rows,  :); late_vals  = late_vals(:);  late_vals  = late_vals(isfinite(late_vals));
[pval, ~, ~] = ranksum(early_vals, late_vals);
fprintf('Threshold-only: p = %.3g (early 0–1 vs late 2–4)\n', pval);

statsCSV = fullfile(fileparts(mfilename('fullpath')), 'stats.csv');
report_stat(statsCSV, 'SFig2F', 'neuron_count_threshold early(0-1) vs late(2-4)', early_vals, late_vals);

% Optional: comparison bar
starStr = 'NS';
if pval < 0.001
    starStr = '***';
elseif pval < 0.01
    starStr = '**';
elseif pval < 0.05
    starStr = '*';
end

y_max = max([curAve(:); mergingTraces(:)], [], 'omitnan');
draw_comparison_bar(y_max * 1.05, y_max * 1.10, {starStr}, blue, 'LineSpacing', 0, 'TextOffset', 0);

% Export
out_w = 4; out_h = 3;
fig2svg("SFig2F", out_w, out_h);
