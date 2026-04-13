%% Build weekDensityData at THRESHOLD from populationDistances

dDiv         = 7;      % days per week
min_neurons  = 5;      % minimum neurons to include a sample
num_shells   = 8;     
shell_width  = 0.1;    % mm (100 µm shells: 0–0.1, 0.1–0.2, ..., 0.7–0.8)
densCells = cell(5, num_shells);   % weeks 0..4 -> 1..5

numbers = [];  % (optional) track neuron counts per sample

for aInd = 1:numel(anID)
    uChan               = uChanAll{aInd};              % channels
    uChanCurrs          = uChanCurrsRef{aInd};         % [Nmaps x 2] = [channel, current]
    daysTrained         = daysTrainedAll{aInd}(:);     % days
    popDates            = popDatesAll{aInd};           % {p} -> vector of day indices
    populationDistances = populationDistancesAll{aInd};% {p} -> {per-day distances}
    daysThresholds      = daysThresholdsAll{aInd};     % [ch x day] threshold in µA

    if isempty(uChan) || isempty(uChanCurrs), continue; end

    % Loop channels
    for chInd = 1:numel(uChan)

        % Loop days
        for dIdx = 1:numel(daysTrained)
            curThresh = daysThresholds(chInd, dIdx);

            % Require valid threshold; optionally constrain range (e.g. 4-6 µA)
            if ~isfinite(curThresh)
                continue;
            end
            % if curThresh < 4 || curThresh > 6
            %     continue;
            % end

            % Find maps that match this channel & this threshold current
            hit = find(uChanCurrs(:,1) == uChan(chInd) & ...
                       uChanCurrs(:,2) == round(curThresh));  % round if mismatch in precision
            if isempty(hit), continue; end

            % For each matching map, see if it has this day
            for pIdx = hit(:).'
                nDates = popDates{pIdx};
                if isempty(nDates), continue; end

                di = find(nDates == dIdx, 1, 'first');
                if isempty(di), continue; end

                dvec_um = populationDistances{pIdx}{di};
                if isempty(dvec_um), continue; end

                % Clean distances
                dvec_um = dvec_um(:);
                dvec_um = dvec_um(isfinite(dvec_um) & dvec_um > 0 & dvec_um <= 800);
                if numel(dvec_um) < min_neurons
                    continue;
                end
                numbers = [numbers; numel(dvec_um)];

                % Week bin (0..4 only)
                wk = ceil(daysTrained(dIdx) / dDiv);
                if ~isfinite(wk) || wk < 0 || wk > 4
                    continue;
                end
                wkSlot = wk + 1;   % 1..5

                % Convert to mm
                dvec_mm = dvec_um / 1000;

                % Compute shell densities
                rho_shell = NaN(num_shells, 1);
                for s = 1:num_shells
                    r_in  = shell_width * (s - 1);
                    r_out = shell_width * s;
                    n_s   = sum(dvec_mm > r_in & dvec_mm <= r_out);
                    V_s   = (4/3) * pi * (r_out^3 - r_in^3);
                    if V_s > 0
                        rho_shell(s) = n_s / V_s;
                    end
                end

                % Append this sample’s shell densities into densCells
                for s = 1:num_shells
                    if isfinite(rho_shell(s))
                        densCells{wkSlot, s} = [densCells{wkSlot, s}; rho_shell(s)];
                    end
                end
            end
        end
    end
end

% ---- Reduce densCells -> median/mean/std per (week, shell)

weekDensityData_median = NaN(5, num_shells);  % weeks x shells
weekDensityData_mean   = NaN(5, num_shells);
weekDensityData_std    = NaN(5, num_shells);

for wk = 1:5
    for s = 1:num_shells
        v = densCells{wk, s};
        if ~isempty(v)
            v = v(isfinite(v));   % keep zeros, just drop NaNs
            if ~isempty(v)
                weekDensityData_median(wk, s) = median(v, 'omitnan');
                weekDensityData_mean(  wk, s) = mean(  v, 'omitnan');
                weekDensityData_std(   wk, s) = std(   v, [], 'omitnan');
            end
        end
    end
end

% Use median as input to the ring plot
weekDensityData = weekDensityData_median;   % [5 x num_shells], week 0..4

%% Ring plot of threshold activation profiles (Fig 2-style)
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
       'PaperUnits','inches','PaperSize',[W H], ...
       'PaperPosition',[0 0 W H]);
tiledlayout(1,5,'TileSpacing','compact','Padding','compact');

cmap = jet(256);
mapIdx = @(v,mx) max(1, min(256, floor(255*(min(v,mx)./mx)) + 1));

for t = 1:5  % weeks 0..4
    nexttile; hold on

    vals = weekDensityData(t, :).';   % [num_shells x 1]
    if all(isnan(vals))
        axis off
        title("Week " + (t-1), 'FontWeight', 'normal');
        continue;
    end

    zeroMask = (vals == 0 | isnan(vals));
    idx      = mapIdx(vals, globalCap);
    colors   = cmap(idx, :);

    for k = num_shells:-1:1
        x1 = 0.5 * abs(k - (num_shells + 1));  % center shift for each ring
        pos = [x1 x1 k k];

        if zeroMask(k)
            fc = cmap(1,:);    % lowest color for empty/no-data shells
        else
            fc = colors(k,:);
        end

        rectangle('Position', pos, 'Curvature', [1 1], ...
                  'FaceColor', fc);
    end
    axis square
    xlim([0 num_shells+1]);
    ylim([0 num_shells+1]);
    axis off
    title("Week " + (t-1), 'FontWeight', 'normal');
end

clim([0 globalCap]); 
colormap(jet);
cb = colorbar; cb.Layout.Tile = 'east';
yl = ylabel(cb,'Neurons / mm3');
yl.Position(1) = min(xlim(cb)) - 1;
cb.Position(3) = 0.1 * cb.Position(3); 
cb.Ticks = [0, 1000, 2000];

sgtitle('At threshold', fontsize=4);
out_w = 6; % inches
out_h = 2; % inches
fig2svg("Fig2F", out_w, out_h);