%% Fig 2c
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
load(FIG2_DATA);

fontsize = 8;
set(groot,'DefaultTextInterpreter','tex',...
          'DefaultAxesTickLabelInterpreter','tex',...
          'DefaultLegendInterpreter','tex',...
          'DefaultTextFontName','Arial',...
          'DefaultAxesFontName','Arial',...
          'DefaultLegendFontName','Arial',...
          'DefaultAxesFontSize', fontsize, ...   
            'DefaultTextFontSize', fontsize, ...  
            'DefaultLegendFontSize', fontsize, ...
            'DefaultColorbarFontSize', fontsize);  

blue   = [0.298, 0.447, 0.690];
orange = [0.867, 0.518, 0.322];
green  = [0.333, 0.659, 0.408];
colors = [blue; orange; green];

set(groot, 'DefaultLegendAutoUpdate', 'off');
figure(); hold on
ax = gca;

uCurrOverall = unique(cell2mat(uCurrAll));
uCurrOverall = uCurrOverall(uCurrOverall>=4 & uCurrOverall<=6);
uChanOverall = unique(cell2mat(uChanAll));
daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
daysTrainedOverall = unique(ceil(daysTrainedOverall/7));

global_b_weeks = [];
global_b_counts = [];

for cInd = 1:numel(uCurrOverall) % process each current

    % Structures for averaging current specific response
    mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
    traceNum=0;

    % Process each animal
    for aInd = 1:numel(anID)
        % if aInd ~= 6; continue; end
        uChan = uChanAll{aInd};
        uCurr = uCurrAll{aInd};
        groupedCurrChans = CurrChansAll{aInd};
        daysTrained = daysTrainedAll{aInd};
        uChanCurrs = uChanCurrsRef{aInd};
        popDates   = popDatesAll{aInd};
        populationDistances = populationDistancesAll{aInd};

        curInd = find(uCurr==uCurrOverall(cInd));
        if isempty(curInd), continue; end
        for chInd = 1:numel(uChan) % process each channel
            if(aInd==5 && chInd==3); continue; end % Skip ICMS 100 channel 6
            neuCnt = NaN(numel(daysTrained), 1);
            hit = find(uChanCurrs(:,1) == uChan(chInd) & uChanCurrs(:,2) == uCurrOverall(cInd));
            if ~isempty(hit)
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
            end
            neuCnt = neuCnt(:);
            valid  = ~isnan(neuCnt) & neuCnt~=0;
            neuCnt  = neuCnt(valid);
            curWeeks = ceil(daysTrained(valid)/7);
            if isempty(neuCnt); continue; end
            % curWeeks = ceil(curDays/dDiv);
            % Group values from same week together
            uW = unique(curWeeks);
            tNeuCnt = zeros(numel(uW),1);
            for w = 1:numel(uW)
               tNeuCnt(w) = median(neuCnt(curWeeks==uW(w)));
            end
            neuCnt = tNeuCnt;
            curDays = uW;  
            scatter(curDays,neuCnt,15,colors(cInd, :),'filled', 'MarkerFaceAlpha', 0.7);
            global_b_weeks = [global_b_weeks; curDays(:)];
            global_b_counts = [global_b_counts; neuCnt(:)];
            % Add trace to current average
            traceNum = traceNum+1;
            alignedDays = zeros(numel(curDays),1);
            for d = 1:numel(curDays)
                alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
            end
            mergingTraces(alignedDays,traceNum) = neuCnt;
        end
    end
    % Average response for current and plot
    curAve = median(mergingTraces,2,'omitnan');
    weekGrid = daysTrainedOverall; 
    plotMask = (weekGrid >= 0 & weekGrid <= 4);
    weeks_to_plot = weekGrid(plotMask);
    num_neurons = curAve(plotMask);
    xticks(0:4);

    early_rows = (weekGrid >= 0 & weekGrid <= 1);
    late_rows  = (weekGrid >= 2 & weekGrid <= 4);
    
    early_vals = mergingTraces(early_rows, :); early_vals = early_vals(:);
    late_vals  = mergingTraces(late_rows,  :); late_vals  = late_vals(:);

    early_vals = early_vals(~isnan(early_vals));
    late_vals  = late_vals(~isnan(late_vals));
    
    [pval, h, stats] = ranksum(early_vals, late_vals);
    fprintf('Current %d µA: p = %.3g\n', uCurrOverall(cInd), pval);

    statsCSV = fullfile(fileparts(mfilename('fullpath')), 'stats.csv');
    report_stat(statsCSV, 'Fig2C', sprintf('neuron_count_%duA early(0-1) vs late(2-4)', uCurrOverall(cInd)), early_vals, late_vals);
  
    p1 = plot(weeks_to_plot, num_neurons);
    tcolor = colors(cInd, :);
    p1.Color = tcolor(1:3);
    p1.LineWidth=2;
end

ylabel('Number of neurons')
xlabel('Weeks of training')
hold on
p4 = patch([-1 -1 -1 -1], [0 0 1 1], blue, 'EdgeColor', 'none');
p5 = patch([-1 -1 -1 -1], [0 0 1 1], orange, 'EdgeColor', 'none'); 
p6 = patch([-1 -1 -1 -1], [0 0 1 1], green, 'EdgeColor', 'none'); 

legend([p4 p5 p6], {'4 µA','5 µA','6 µA'}, 'Location', 'northwest', 'Interpreter', 'tex', 'Box', 'off');
hold off
xlim([-0.5, 4.5]); ylim([-10, 1220])
draw_comparison_bar(920, 970, {'*', '*', '**'}, colors, 'fontsize', 20, 'LineSpacing', 50, 'TextOffset', -60);
out_w = 4;          % inches
out_h = 3;          % inches
fig2svg("Fig2C", out_w, out_h);
% save('\\10.129.151.108\xieluanlabs\xl_stimulation\Roy\robin - Copy\figures\SFig2D\behavior_count_source.mat', 'global_b_weeks', 'global_b_counts');

