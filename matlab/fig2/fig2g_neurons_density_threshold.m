%% Fig 2G - Number of neurons (left) and activation density (right) at threshold
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

dDiv = 7;  % days per week bin

% --- Collectors ---
weeks           = [];
allNeuronCount  = [];   % number of neurons (within 800 µm)
allDensity0_200 = [];   % density within 0-200 µm sphere (neurons/mm^3)

for aInd = 1:numel(anID)
    uChan               = uChanAll{aInd};
    daysTrained         = daysTrainedAll{aInd};
    weeksTrained        = ceil(daysTrained./7);
    uChanCurrs          = uChanCurrsRef{aInd};
    popDates            = popDatesAll{aInd};
    populationDistances = populationDistancesAll{aInd};
    daysThresholds      = daysThresholdsAll{aInd};

    for chInd = 1:numel(uChan)
        for cDay = 1:numel(daysTrained)
            curThresh = round(daysThresholds(chInd, cDay));
            if ~isfinite(curThresh), continue; end

            hit = find(uChanCurrs(:,1)==uChan(chInd) & uChanCurrs(:,2)==curThresh);
            if isempty(hit), continue; end

            for pIdx = hit(:).'
                nDates = popDates{pIdx};
                dInd   = find(nDates==cDay, 1, 'first');
                if isempty(dInd), continue; end

                dvec_um = populationDistances{pIdx}{dInd};
                if isempty(dvec_um), continue; end
                dvec_um = dvec_um(:);
                dvec_um = dvec_um(isfinite(dvec_um) & dvec_um > 0 & dvec_um <= 800);
                if numel(dvec_um) < 1, continue; end

                % Neuron count
                nCount = numel(dvec_um);

                % Density within 0-200 µm
                dvec_mm    = dvec_um / 1000;
                r_fixed_mm = 0.2;
                n_0_200    = sum(dvec_mm <= r_fixed_mm);
                vol_0_200  = (4/3)*pi*(r_fixed_mm^3);
                rho_0_200  = n_0_200 / vol_0_200;

                weeks           = [weeks;           weeksTrained(cDay)];
                allNeuronCount  = [allNeuronCount;  nCount];
                allDensity0_200 = [allDensity0_200; rho_0_200];
            end
        end
    end
end

% Filter
mask = isfinite(weeks) & weeks>=0 & weeks<=4 & ...
       isfinite(allNeuronCount) & allNeuronCount>0 & ...
       isfinite(allDensity0_200) & allDensity0_200>0;

weeks           = weeks(mask);
allNeuronCount  = allNeuronCount(mask);
allDensity0_200 = allDensity0_200(mask);

% Per-week medians
medCount   = accumarray(weeks+1, allNeuronCount,  [5 1], @(x) median(x,'omitnan'), NaN);
medDensity = accumarray(weeks+1, allDensity0_200, [5 1], @(x) median(x,'omitnan'), NaN);

wk_has_cnt = ~isnan(medCount);
wk_has_den = ~isnan(medDensity);
x_cnt = find(wk_has_cnt)-1 - 0.1;
x_den = find(wk_has_den)-1 + 0.1;

% --- Colors ---
col_count   = [0.506, 0.447, 0.702];  % seaborn deep purple = neuron count (left)
col_density = [0.839, 0.153, 0.157];  % seaborn deep red = density (right)

% --- Plot ---
f = figure('Position',[2000 200 800 500],'Renderer','painters'); hold on
set(findall(f,'-property','Interpreter'),'Interpreter','none');
set(findall(f,'-property','FontName'),'FontName','Arial');

% LEFT axis - neuron count
yyaxis left; axL = gca;
scatter(weeks - 0.1, allNeuronCount, 18, 'filled', ...
    'MarkerFaceColor', col_count, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6);
plot(x_cnt, medCount(wk_has_cnt), '-', 'Color', col_count, 'LineWidth', 1.5);
ylabel('Number of neurons');
axL.YColor = col_count;

% RIGHT axis - density
yyaxis right; axR = gca;
scatter(weeks + 0.1, allDensity0_200, 18, 'filled', ...
    'MarkerFaceColor', col_density, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6);
plot(x_den, medDensity(wk_has_den), '-', 'Color', col_density, 'LineWidth', 1.5);
ylabel('Activation density (neurons / mm^3)', 'Rotation', 270);
axR.YColor = col_density;

xlabel('Weeks of training'); xlim([-0.5 4.5]); xticks(0:4); box off; grid off

% Set Y limits
yyaxis left;  axL.YLim = [-max(allNeuronCount)*0.05, max(allNeuronCount)*1.2];
yyaxis right; axR.YLim = [0, 6000];

% --- Stats ---
p_count   = ranksum(allNeuronCount(weeks<=1),  allNeuronCount(weeks>=2));
p_density = ranksum(allDensity0_200(weeks<=1), allDensity0_200(weeks>=2));
fprintf('Threshold - p_count=%.3g, p_density=%.3g\n', p_count, p_density);

statsCSV = fullfile(fileparts(mfilename('fullpath')), 'stats.csv');
report_stat(statsCSV, 'Fig2G', 'neuron_count_threshold early(0-1) vs late(2-4)', ...
    allNeuronCount(weeks<=1), allNeuronCount(weeks>=2));
report_stat(statsCSV, 'Fig2G', 'density_threshold early(0-1) vs late(2-4)', ...
    allDensity0_200(weeks<=1), allDensity0_200(weeks>=2));

% Comparison bar
yyaxis right
draw_comparison_bar(axR.YLim(2)*0.85, axR.YLim(2)*0.92, {'NS', 'NS'}, [col_count; col_density], ...
    'LineSpacing', axR.YLim(2)*0.08, 'TextOffset', 0);

fig2svg("Fig2G", 3, 2);
