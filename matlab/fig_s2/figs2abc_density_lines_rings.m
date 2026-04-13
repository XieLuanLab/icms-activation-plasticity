% SFig 2A (activation density vs distance from stim site at 4 uA)

% Build weekDensityData from populationDistances
dDiv = 7;  % days per week

targetCurrents = [4,5]; 
C = numel(targetCurrents);
min_neurons = 1;
num_shells = 20;    
shell_width = 0.05;   % mm (50 µm)
num_shells = 8;    
shell_width = 0.1;   % mm (100 µm)
% Collectors: one cell per (current, week, shell), will hold a vector of samples
densCells = cell(C, 5, num_shells);

numbers = [];
for aInd = 1:numel(anID)
    uChan               = uChanAll{aInd};             % [channels]
    uCurr               = uCurrAll{aInd};             % list of currents in this animal
    daysTrained         = daysTrainedAll{aInd}(:);    % [days]
    uChanCurrs          = uChanCurrsRef{aInd};        % [Nmaps x 2] = [channel, current]
    popDates            = popDatesAll{aInd};          % cell of date indices for each map
    populationDistances = populationDistancesAll{aInd}; % cell of { per-day vector of distances in µm }

    if isempty(uChan) || isempty(uCurr) || isempty(uChanCurrs), continue; end

    % For each desired current
    for cIdx = 1:C
        curCurr = targetCurrents(cIdx);

        % quick skip if this animal doesn’t have this current
        if ~ismember(curCurr, uCurr), continue; end

        % Loop channels
        for chInd = 1:numel(uChan)
            % Find all "maps" (rows) for this channel & this current
            hit = find(uChanCurrs(:,1) == uChan(chInd) & uChanCurrs(:,2) == curCurr);
            if isempty(hit), continue; end

            % Each hit corresponds to a mapping table with its own dates
            for pIdx = hit(:).'
                nDates = popDates{pIdx};  % vector of day indices into daysTrained
                if isempty(nDates), continue; end

                % Loop days within this mapping
                for di = 1:numel(nDates)
                    dIdx = nDates(di);                     % 1-based index into daysTrained
                    if ~isfinite(dIdx) || dIdx < 1 || dIdx > numel(daysTrained), continue; end

                    wk = ceil(daysTrained(dIdx) / dDiv);   % week bin
                    if ~isfinite(wk) || wk < 0 || wk > 4, continue; end

                    dvec_um = populationDistances{pIdx}{di};  % distances (µm) for this day
                    if isempty(dvec_um), continue; end

                    % Clean and require at least 1 distance to form shells
                    dvec_um = dvec_um(:);
                    dvec_um = dvec_um(isfinite(dvec_um) & dvec_um > 0 & dvec_um <= 800);
                    if numel(dvec_um) < min_neurons, continue; end
                    number = numel(dvec_um);
                    numbers = [numbers; number];

                    % Convert to mm
                    dvec_mm = dvec_um / 1000;

                    % Compute shell densities for 10 shells of 0.1 mm
                    rho_shell = NaN(num_shells, 1);
                    for s_index = 1:num_shells
                        r_in  = shell_width * (s_index - 1);
                        r_out = shell_width * s_index;
                        n_s   = sum(dvec_mm > r_in & dvec_mm <= r_out);
                        V_s   = (4/3) * pi * (r_out^3 - r_in^3);
                        if V_s > 0
                            rho_shell(s_index) = n_s / V_s;
                        end
                    end


                    % Append this sample’s shell densities to the collectors
                    wkSlot = wk + 1;  % store weeks 0..4 at indices 1..5
                    for s = 1:num_shells
                        densCells{cIdx, wkSlot, s} = [densCells{cIdx, wkSlot, s}; rho_shell(s)];
                    end
                end
            end
        end
    end
end

% ---- Reduce cells → median/mean/std arrays
weekDensityData_median = NaN(C, 5, num_shells);
weekDensityData_mean   = NaN(C, 5, num_shells);
weekDensityData_std    = NaN(C, 5, num_shells);

for cIdx = 1:C
    for wk = 1:5
        for s = 1:num_shells
            % v is cell array with each cell containing array of densities 
            v = densCells{cIdx, wk, s};
            if ~isempty(v)
                % keep zeros (don’t filter them out) so it matches ring-plot semantics
                vv = v(isfinite(v));
                if ~isempty(vv)
                    % summary statistic of cell array stored in 3d array 
                    weekDensityData_median(cIdx, wk, s) = median(vv, 'omitnan');
                    weekDensityData_mean(  cIdx, wk, s) = mean(  vv, 'omitnan');
                    weekDensityData_std(   cIdx, wk, s) = std(   vv, [], 'omitnan');
                end
            end
        end
    end
end

% pick which to feed the ring plot
weekDensityData = weekDensityData_median;


%% 
fontsize = 10;
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

% x-values for shells (midpoints, µm)
shell_edges_um = 0:100:800;     % 8 shells (0–800 µm)
x = shell_edges_um(2:end) - 50;  % midpoints: 50,150,...,950

% seaborn deep color palette
deep_colors = [
    0.298, 0.447, 0.690;   % blue
    0.866, 0.517, 0.321;   % orange
    0.333, 0.659, 0.408;   % green
    0.769, 0.306, 0.322;   % red
    0.506, 0.447, 0.702;   % purple
];
%% SFig2A
figure; hold on
c = find(targetCurrents == 4);

for t = 1:5
    vals = squeeze(weekDensityData(c,t,:));  % [10×1]
    plot(x, vals, 'Color', deep_colors(t,:), 'LineWidth', 1.0);  % thinner line
end

xlabel('Distance from stimulation site (μm)');
ylabel('Activation density (neurons / mm^3)');

legend({'Week 0','Week 1','Week 2','Week 3','Week 4'}, ...
    'Box','off', ...
    'Location','northeast', ...
    'FontSize', fontsize - 1);

xlim([0 800]);
ylim([0 inf]);
xticks(0:200:800);             % fixed 200 µm tick labels
xticklabels({'0','200','400','600','800'}); % stable across export

box off;
set(gca,'TickDir','out','LineWidth',0.75);
title('At 4 μA')

out_w = 4; % inches
out_h = 2; % inches
fig2svg("SFig_4uA_density_line", out_w, out_h);

%% SFig2B
figure; hold on
c = find(targetCurrents == 5);

for t = 1:5
    vals = squeeze(weekDensityData(c,t,:)); 
    plot(x, vals, 'Color', deep_colors(t,:), 'LineWidth', 1.0);  % thinner line
end

xlabel('Distance from stimulation site (μm)');
ylabel('Activation density (neurons / mm^3)');

legend({'Week 0','Week 1','Week 2','Week 3','Week 4'}, ...
    'Box','off', ...
    'Location','northeast', ...
    'FontSize', fontsize - 1);

xlim([0 800]);
ylim([0 inf]);
xticks(0:200:800);             % fixed 200 µm tick labels
xticklabels({'0','200','400','600','800'}); % stable across export

box off;
set(gca,'TickDir','out','LineWidth',0.75);
title('At 5 μA')

out_w = 4; % inches
out_h = 2; % inches
fig2svg("SFig_5uA_density_line", out_w, out_h);

%% Fig 2C
fontsize = 5;
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

ring_in = 5.0;  pad_in = 0.10;
W = 5*ring_in + 2*pad_in;  
H = 1*ring_in + 2*pad_in; 

figure('Units','inches','Position',[22 5 W H], ...
       'PaperUnits','inches','PaperSize',[W H],'PaperPosition',[0 0 W H]);
tiledlayout(1,5,'TileSpacing','compact','Padding','compact');

cmap = jet(256);
current_uA = 4;
c = find(targetCurrents == current_uA); 
fprintf('Current set to %d uA\n', current_uA);

%  Global cap - match Fig2D colorbar (hard-cap at 2000)
if ~exist('globalCap', 'var')
    allVals   = weekDensityData;
    finiteVals = allVals(isfinite(allVals));
    if isempty(finiteVals), error('weekDensityData is empty/NaN'); end
    globalCap = min(prctile(finiteVals, 99), 2000);
end
globalCap = min(globalCap, 2000);  % enforce 2000 ceiling
% Helper: saturate and map to 1..256 (v=cap -> 256)
mapIdx = @(v,mx) max(1, min(256, floor(255*(min(v,mx)./mx)) + 1));

for t = 1:5
    nexttile; hold on
    vals = squeeze(weekDensityData(c,t,:));   
    if ~isnan(vals(1))
        zeroMask = (vals == 0);
        idx = mapIdx(vals, globalCap);
        colors = cmap(idx, :);

        for k = num_shells:-1:1
            % Compute a visually consistent radius offset
            % x1 ~ center shift for each ring
            x1 = 0.5 * abs(k - (num_shells + 1));
            pos = [x1 x1 k k];

            if zeroMask(k)
                fc = cmap(1,:);   % lowest color for empty shells
            else
                fc = colors(k,:);
            end
            rectangle('Position', pos, 'Curvature', [1 1], 'FaceColor', fc)
        end
    end

    axis square
    xlim([0 num_shells+1])
    ylim([0 num_shells+1])
    axis off
    title("Week " + (t-1), 'FontWeight', 'normal');
end

clim([0 2000]);  % fixed to match Fig2D
colormap(jet);
cb = colorbar; cb.Layout.Tile = 'east';
yl = ylabel(cb,'Neurons / mm3');
yl.Position(1) = min(xlim(cb)) - 1;
cb.Position(3) = 0.1 * cb.Position(3); 
cb.Ticks = [0, 1000, 2000];
out_w = 6; % inches
out_h = 2; % inches
fig2svg("SFig_4uA_rings", out_w, out_h);
