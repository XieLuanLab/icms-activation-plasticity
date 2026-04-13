%% SFig2D - Activation density: Control vs Behavioral at 5 µA (800 µm truncation)
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
load(SFIG2D_CONTROL);

Xd = gR;
Yd = allR200;
valid = ~isnan(Xd) & ~isnan(Yd) & Xd < 4 & Yd ~= 0;
Xd = Xd(valid);  Yd = Yd(valid);

countData = struct([]);
xu_den = unique(Xd);
for k = 1:numel(xu_den)
    idx = (Xd == xu_den(k));
    countData(k).week = xu_den(k);
    countData(k).x    = Xd(idx);
    countData(k).y    = Yd(idx);
end

%% Generate behavioral density at 5 µA with 800 µm truncation + full Fig2E mask, weeks 0-3
curCurr = 5; min_neurons = 1;
wks=[]; allDistUm=[]; allRhoR50=[]; allRho200=[];
allShellMean=[]; allShellMedian=[];

for aInd = 1:numel(anID)-1
    uChan=uChanAll{aInd}; uCurr=uCurrAll{aInd};
    daysTrained=daysTrainedAll{aInd};
    uChanCurrs=uChanCurrsRef{aInd};
    popDates=popDatesAll{aInd};
    populationDistances=populationDistancesAll{aInd};
    curInd=find(uCurr==curCurr,1);
    if isempty(curInd), continue; end
    for chInd=1:numel(uChan)
        hit=find(uChanCurrs(:,1)==uChan(chInd)&uChanCurrs(:,2)==curCurr);
        if isempty(hit), continue; end
        for pIdx=hit(:).'
            nDates=popDates{pIdx}; curDays=daysTrained(nDates);
            weeksBin=ceil(curDays/7);
            for di=1:numel(weeksBin)
                wk=weeksBin(di);
                if ~isfinite(wk)||wk<0||wk>3, continue; end
                dvec_um=populationDistances{pIdx}{di};
                if isempty(dvec_um), continue; end
                dvec_um=dvec_um(:);
                dvec_um=dvec_um(isfinite(dvec_um)&dvec_um>0&dvec_um<=800);
                if numel(dvec_um)<min_neurons, continue; end
                r50_um=median(dvec_um,'omitnan');
                r50_mm=r50_um/1000; vol_mm3=(4/3)*pi*(r50_mm^3);
                if vol_mm3==0, continue; end
                rho_r50=sum(dvec_um<=r50_um)/vol_mm3;
                dvec_mm=dvec_um/1000;
                rho_200=sum(dvec_mm<=0.2)/((4/3)*pi*(0.2^3));
                rho_shell=NaN(5,1);
                for s=1:5
                    r_in=0.1*(s-1); r_out=0.1*s;
                    n_s=sum(dvec_mm>r_in&dvec_mm<=r_out);
                    V_s=(4/3)*pi*(r_out^3-r_in^3);
                    rho_shell(s)=n_s/V_s;
                end
                wks=[wks;wk]; allDistUm=[allDistUm;r50_um];
                allRhoR50=[allRhoR50;rho_r50]; allRho200=[allRho200;rho_200];
                allShellMean=[allShellMean;mean(rho_shell,'omitnan')];
                allShellMedian=[allShellMedian;median(rho_shell,'omitnan')];
            end
        end
    end
end

% Full Fig2E mask
mask = isfinite(wks)&wks>=0&wks<=3& ...
       isfinite(allDistUm)&allDistUm>0& ...
       isfinite(allRhoR50)&allRhoR50>0& ...
       isfinite(allRho200)&allRho200>0& ...
       isfinite(allShellMean)&allShellMean>0& ...
       isfinite(allShellMedian)&allShellMedian>0;

b_allWeeks = wks(mask);
b_allVals  = allRho200(mask);

%% Build behavioral struct
b_weeks_unique = unique(b_allWeeks);
b_densityData = struct([]);
for k = 1:numel(b_weeks_unique)
    idx = (b_allWeeks == b_weeks_unique(k));
    b_densityData(k).week = b_weeks_unique(k);
    b_densityData(k).x    = b_allWeeks(idx);
    b_densityData(k).y    = b_allVals(idx);
end

nWeeks = min(numel(countData), numel(b_densityData));

%% Plot
col_control  = [0.333, 0.659, 0.408];
col_behavior = [0.506, 0.447, 0.702];
colors = [col_control; col_behavior];

figure; hold on

allWeeks = []; allVals = [];
for k = 1:nWeeks
    wk = countData(k).week; ys = countData(k).y;
    allWeeks = [allWeeks; wk*ones(numel(ys),1)];
    allVals  = [allVals; ys(:)];
end

b_allWeeks_plot = []; b_allVals_plot = [];
for k = 1:nWeeks
    wk = b_densityData(k).week; ys = b_densityData(k).y;
    b_allWeeks_plot = [b_allWeeks_plot; wk*ones(numel(ys),1)];
    b_allVals_plot  = [b_allVals_plot; ys(:)];
end

s1 = scatter(allWeeks - 0.1, allVals, 24, 'filled', ...
    'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor',col_control);
s2 = scatter(b_allWeeks_plot + 0.1, b_allVals_plot, 24, 'filled', ...
    'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6, 'MarkerFaceColor',col_behavior);

% Legend - create immediately after scatter, before anything else
legend([s1, s2], {'Control', 'Behavior'}, ...
    'Location','northwest', 'Box','off', 'AutoUpdate', 'off');

weeks_unique = unique(allWeeks);
weekMedian = arrayfun(@(wk) median(allVals(allWeeks==wk),'omitnan'), weeks_unique);
b_wu = unique(b_allWeeks_plot);
b_weekMedian = arrayfun(@(wk) median(b_allVals_plot(b_allWeeks_plot==wk),'omitnan'), b_wu);

h1 = plot(weeks_unique - 0.1, weekMedian, 'Color', col_control, 'LineWidth', 1.5);
h2 = plot(b_wu + 0.1, b_weekMedian, 'Color', col_behavior, 'LineWidth', 1.5);
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

xlabel('Weeks of training');
ylabel('Activation density (neurons / mm^3)');
xticks(weeks_unique);
box off
xlim([-0.75, 3.75]);
ylim([-100, 7000])

%% Stats
p_count = ranksum(allVals(allWeeks<=1), allVals(allWeeks>=2));
fprintf('control p_density=%.3g\n', p_count);
b_p_count = ranksum(b_allVals_plot(b_allWeeks_plot<=1), b_allVals_plot(b_allWeeks_plot>=2));
fprintf('behavior p_density=%.3g\n', b_p_count);

statsCSV = fullfile(fileparts(mfilename('fullpath')), 'stats.csv');
report_stat(statsCSV, 'SFig2D', 'control_density early(0-1) vs late(2-3)', ...
    allVals(allWeeks<=1), allVals(allWeeks>=2));
report_stat(statsCSV, 'SFig2D', 'behavior_density_5uA early(0-1) vs late(2-3)', ...
    b_allVals_plot(b_allWeeks_plot<=1), b_allVals_plot(b_allWeeks_plot>=2));

draw_comparison_bar(5500, 6000, {'NS', '*'}, colors, 'fontsize', 12, ...
    'LineSpacing', 200, 'TextOffset', 0, 'LateX', [2, 3]);

out_w = 4; out_h = 3;
fig2svg("SFig2D", out_w, out_h);
