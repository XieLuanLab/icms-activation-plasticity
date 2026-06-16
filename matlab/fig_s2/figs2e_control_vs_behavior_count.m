%% SFig2E - Number of neurons: Control vs Behavioral at 5 µA (800 µm truncation)
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

addpath(fullfile(pwd, 'matlab'));

%% Load control data
load(SFIG2E_CONTROL);

Xd = g;
Yd = allCounts;
valid = ~isnan(Xd) & ~isnan(Yd) & Xd < 5 & Yd ~= 0;
Xd = Xd(valid);  Yd = Yd(valid);

countData = struct([]);
xu_den = unique(Xd);
for k = 1:numel(xu_den)
    idx = (Xd == xu_den(k));
    countData(k).week   = xu_den(k);
    countData(k).x      = Xd(idx);
    countData(k).y      = Yd(idx);
end

%% Generate behavioral count at 5 µA with 800 µm truncation, weeks 0-3
curCurr = 5;
b_allWeeks = [];
b_allVals = [];

for aInd = 1:numel(anID)-1
    uChan = uChanAll{aInd};
    uCurr = uCurrAll{aInd};
    daysTrained = daysTrainedAll{aInd};
    uChanCurrs = uChanCurrsRef{aInd};
    popDates = popDatesAll{aInd};
    populationDistances = populationDistancesAll{aInd};

    curInd = find(uCurr == curCurr);
    if isempty(curInd), continue; end

    for chInd = 1:numel(uChan)
        hit = find(uChanCurrs(:,1) == uChan(chInd) & uChanCurrs(:,2) == curCurr);
        if isempty(hit), continue; end

        neuCnt = NaN(numel(daysTrained), 1);
        for pIdx = hit(:).'
            nDates = popDates{pIdx};
            for di = 1:numel(nDates)
                dIdx = nDates(di);
                dvec_um = populationDistances{pIdx}{di};
                if isempty(dvec_um), continue; end
                dvec_um = dvec_um(:);
                dvec_um = dvec_um(isfinite(dvec_um) & dvec_um > 0 & dvec_um <= 800);
                neuCnt(dIdx) = numel(dvec_um);
            end
        end

        valid_nc = ~isnan(neuCnt) & neuCnt ~= 0;
        neuCnt = neuCnt(valid_nc);
        curWeeks = ceil(daysTrained(valid_nc) / 7);

        % Keep weeks 0-3 only
        keep = curWeeks <= 3;
        neuCnt = neuCnt(keep);
        curWeeks = curWeeks(keep);
        if isempty(neuCnt), continue; end

        % Per-week median for this channel
        uW = unique(curWeeks);
        for w = 1:numel(uW)
            b_allWeeks = [b_allWeeks; uW(w)];
            b_allVals = [b_allVals; median(neuCnt(curWeeks == uW(w)))];
        end
    end
end

%% Build behavioral struct to match control format
b_weeks_unique = unique(b_allWeeks);
b_countData = struct([]);
for k = 1:numel(b_weeks_unique)
    idx = (b_allWeeks == b_weeks_unique(k));
    b_countData(k).week = b_weeks_unique(k);
    b_countData(k).x    = b_allWeeks(idx);
    b_countData(k).y    = b_allVals(idx);
end

% Trim control to same number of weeks as behavioral
nWeeks = min(numel(countData), numel(b_countData));

%% Plot
col_control  = [0.333, 0.659, 0.408];
col_behavior = [0.506, 0.447, 0.702];
colors = [col_control; col_behavior];

figure; hold on

% Collect all values for control
allWeeks = [];
allVals  = [];
for k = 1:nWeeks
    wk = countData(k).week;
    ys = countData(k).y;
    allWeeks = [allWeeks; wk*ones(numel(ys),1)];
    allVals  = [allVals; ys(:)];
end

% Collect all values for behavioral
b_allWeeks_plot = [];
b_allVals_plot  = [];
for k = 1:nWeeks
    wk = b_countData(k).week;
    ys = b_countData(k).y;
    b_allWeeks_plot = [b_allWeeks_plot; wk*ones(numel(ys),1)];
    b_allVals_plot  = [b_allVals_plot; ys(:)];
end

% Scatter
s1 = scatter(allWeeks - 0.1, allVals, 24, 'filled', ...
    'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor',col_control);
s2 = scatter(b_allWeeks_plot + 0.1, b_allVals_plot, 24, 'filled', ...
    'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor',col_behavior);

% Legend - create immediately after scatter, before anything else
legend([s1, s2], {'Control', 'Behavior'}, ...
    'Location','northwest', 'Box','off', 'AutoUpdate', 'off');

% Median lines
weeks_unique = unique(allWeeks);
weekMedian = arrayfun(@(wk) median(allVals(allWeeks==wk), 'omitnan'), weeks_unique);
b_wu = unique(b_allWeeks_plot);
b_weekMedian = arrayfun(@(wk) median(b_allVals_plot(b_allWeeks_plot==wk), 'omitnan'), b_wu);

h1 = plot(weeks_unique - 0.1, weekMedian, 'Color', col_control, 'LineWidth', 1.5);
h2 = plot(b_wu + 0.1, b_weekMedian, 'Color', col_behavior, 'LineWidth', 1.5);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

xlabel('Weeks of training');
ylabel('Number of neurons');
xticks(weeks_unique);
box off
xlim([-0.75, 3.75]);
ylim([-20, 900])

%% Stats
p_count = ranksum(allVals(allWeeks<=1), allVals(allWeeks>=2));
fprintf('control p_count=%.3g\n', p_count);
b_p_count = ranksum(b_allVals_plot(b_allWeeks_plot<=1), b_allVals_plot(b_allWeeks_plot>=2));
fprintf('behavior p_count=%.3g\n', b_p_count);

statsCSV = fullfile(fileparts(mfilename('fullpath')), 'stats.csv');
report_stat(statsCSV, 'SFig2E', 'control_count early(0-1) vs late(2-3)', ...
    allVals(allWeeks<=1), allVals(allWeeks>=2));
report_stat(statsCSV, 'SFig2E', 'behavior_count_5uA early(0-1) vs late(2-3)', ...
    b_allVals_plot(b_allWeeks_plot<=1), b_allVals_plot(b_allWeeks_plot>=2));

% Comparison bar
draw_comparison_bar(700, 750, {'NS', '**'}, colors, 'fontsize', 12, ...
    'LineSpacing', 50, 'TextOffset', 0, 'LateX', [2, 3]);

out_w = 4; out_h = 3;
fig2svg("SFig2E", out_w, out_h);
