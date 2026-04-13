%% Radius that encloses 50%
% Activation density
% Count # of activation radii which is # of neurons and
% then divide by sphere of radius r50
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
load(FIG2_DATA);

dDiv = 7;
curCurr = 5;
min_neurons = 1;

% collectors (one row per observation: site × day)
weeks                 = [];
allDistancesUm        = [];   % r50 in µm
allDensityR50_simple  = [];   % density within r50 (neurons/mm^3)
allDensity0_200       = [];   % density within 0–500 µm sphere (neurons/mm^3)
allDensityR500_shellMean = []; % mean of shell densities (0–100,...,400–500)
allDensityR500_shellMedian = []; % median of shell densities

for aInd = 1:numel(anID)-1
    uChan              = uChanAll{aInd};
    uCurr              = uCurrAll{aInd};
    daysTrained        = daysTrainedAll{aInd};
    uChanCurrs         = uChanCurrsRef{aInd};          
    popDates           = popDatesAll{aInd};           
    populationDistances= populationDistancesAll{aInd};

    curInd = find(uCurr == curCurr, 1);
    if isempty(curInd), continue; end                   
    for chInd = 1:numel(uChan)
        % mapping rows for this channel+current
        hit = find(uChanCurrs(:,1) == uChan(chInd) & uChanCurrs(:,2) == curCurr);
        if isempty(hit), continue; end

        for pIdx = hit(:).'
            nDates = popDates{pIdx};                    % vector of day indices
            curDays = daysTrained(nDates);
            weeksBin = ceil(curDays / dDiv);            % 0..∞

            % loop days in this mapping
            for di = 1:numel(weeksBin)
                wk = weeksBin(di);

                if ~isfinite(wk) || wk < 0 || wk > 4, continue; end

                dvec_um = populationDistances{pIdx}{di};   % distances (µm)
                if isempty(dvec_um), continue; end
                dvec_um = dvec_um(:);
                dvec_um = dvec_um(isfinite(dvec_um) & dvec_um > 0 & dvec_um <= 800);
                if numel(dvec_um) < min_neurons, continue; end
                r50_um = median(dvec_um, 'omitnan');              

                % neurons within r50 (count)
                n_r50 = sum(dvec_um <= r50_um);

                % simple: density within r50 (neurons/mm^3)
                r50_mm = r50_um / 1000;
                vol_mm3 = (4/3)*pi*(r50_mm^3);
                if vol_mm3 == 0, continue; end
                rho_r50_simple = n_r50 / vol_mm3;

                % simple: density within fixed 0–200 µm 
                dvec_mm   = dvec_um / 1000;
                r_fixed_mm = 0.2;
                n_0_200   = sum(dvec_mm <= r_fixed_mm);
                vol_0_200 = (4/3)*pi*(r_fixed_mm^3);
                rho_0_200 = n_0_200 / vol_0_200;
               
                
                % shell: 5 shells of 0.1 mm each, densities per shell
                % this is a bit misleading since outer shells are larger in
                % volume so will drag down mean/median density 
                rho_shell = NaN(5,1);
                for s = 1:5
                    r_in  = 0.1*(s-1);
                    r_out = 0.1*s;
                    n_s   = sum(dvec_mm > r_in & dvec_mm <= r_out);
                    V_s   = (4/3)*pi*(r_out^3 - r_in^3);
                    rho_shell(s) = n_s / V_s;
                end
                shell_mean   = mean(rho_shell,   'omitnan');
                shell_median = median(rho_shell, 'omitnan');
                
                weeks = [weeks; wk];
                allDistancesUm              = [allDistancesUm; r50_um];
                allDensityR50_simple        = [allDensityR50_simple; rho_r50_simple];
                allDensity0_200             = [allDensity0_200; rho_0_200];
                allDensityR500_shellMean    = [allDensityR500_shellMean;   shell_mean];
                allDensityR500_shellMedian  = [allDensityR500_shellMedian; shell_median];
            end
        end
    end
end

mask = isfinite(weeks) & weeks>=0 & weeks<=4 & ...
       isfinite(allDistancesUm) & allDistancesUm>0 & ...
       isfinite(allDensityR50_simple) & allDensityR50_simple>0 & ...
       isfinite(allDensity0_200) & allDensity0_200>0 & ...
       isfinite(allDensityR500_shellMean) & allDensityR500_shellMean>0 & ...
       isfinite(allDensityR500_shellMedian) & allDensityR500_shellMedian>0;

weeks                      = weeks(mask);
allDistancesUm             = allDistancesUm(mask);
allDensityR50_simple       = allDensityR50_simple(mask);
allDensity0_200            = allDensity0_200(mask);
allDensityR500_shellMean   = allDensityR500_shellMean(mask);
allDensityR500_shellMedian = allDensityR500_shellMedian(mask);

% -------- per-week medians for both metrics
% Choose which density metric to plot
% Options: 'r50', 'sphere500', 'shellMean', 'shellMedian'
method = 'sphere200';


yDensity     = allDensity0_200;
medDensity   = accumarray(weeks+1, allDensity0_200, [], @(x) median(x,'omitnan'));
% --- Plot ---
col_density  = [0.839, 0.153, 0.157];  % seaborn deep red
col_distance = [0, 0, 0];
densityLabel = 'Activation density (neurons / mm^3)';

f = figure('Position',[2000 200 800 500],'Renderer','painters'); hold on
set(findall(f,'-property','Interpreter'),'Interpreter','none');
set(findall(f,'-property','FontName'),'FontName','Arial');
set(findall(f,'-property','FontWeight'),'FontWeight','normal');

yyaxis left
scatter(weeks - 0.1, allDistancesUm, 18, 'filled', ...
    'MarkerFaceColor', col_distance, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6);
plot((0:4) - 0.1, accumarray(weeks+1, allDistancesUm, [], @(x) median(x,'omitnan')), ...
    '-', 'Color', col_distance, 'LineWidth', 1.5);
ylabel('Enclosing radius (µm)');
ax = gca;  ax.YColor = col_distance;

yyaxis right
scatter(weeks + 0.1, yDensity, 18, 'filled', ...
    'MarkerFaceColor', col_density, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6);
plot((0:4) + 0.1, medDensity, '-', 'Color', col_density, 'LineWidth', 1.5);
ylabel(densityLabel, rotation=270);
ax.YColor = col_density;

xlabel('Weeks of training'); 
xlim([-0.5 4.5]); xticks(0:4); box off; grid off;

[p_r50, t, p] = ranksum(allDistancesUm(weeks<=1), allDistancesUm(weeks>=2));
p_den = ranksum(yDensity(weeks<=1), yDensity(weeks>=2 ));
fprintf('Current %.1f µA - p_r50=%.3g, p_density=%.3g [%s]\n', curCurr, p_r50, p_den, method);

statsCSV = fullfile(fileparts(mfilename('fullpath')), 'stats.csv');
report_stat(statsCSV, 'Fig2E', 'r50_5uA early(0-1) vs late(2-4)', allDistancesUm(weeks<=1), allDistancesUm(weeks>=2));
report_stat(statsCSV, 'Fig2E', 'density_5uA early(0-1) vs late(2-4)', yDensity(weeks<=1), yDensity(weeks>=2));

yyaxis right
ylim([0, 6000]);
yyaxis left
ylim([0, 600]);

yyaxis right;
draw_comparison_bar(4900,5100, {'*', '***'}, [col_density; col_distance], 'LineSpacing', 500, 'TextOffset', 0);

ylim([0, 6000])
yyaxis left; 
ylim([0, 600])
fig2svg("Fig2E", 2.5, 2);


for week = 0:4
    y_week = allDistancesUm(weeks == week);
    mu = mean(y_week);
    sigma = std(y_week);
    n = numel(y_week);
    fprintf('Week %d: %.1f +/- %.1f (%d points)\n', week, mu, sigma, n);
end

