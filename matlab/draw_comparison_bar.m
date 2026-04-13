function h = draw_comparison_bar(group_y, y_top, sig_texts, colors, varargin)
% Draws comparison bars with multi-line, colorized labels.
% Always expects colors as N×3 numeric matrix (each row = RGB for a line).
%
% Usage:
%   draw_comparison_bar(0.5, 0.8, {'**','p=0.003','n=40'}, colors);

    % Parse options
    p = inputParser;
    addParameter(p, 'BarColor', [0 0 0]);
    addParameter(p, 'LineWidth', 1.5);
    addParameter(p, 'EarlyX', [0 1]);
    addParameter(p, 'LateX',  [2 4]);
    addParameter(p, 'TextOffset', 0.02);
    addParameter(p, 'LineSpacing', 0.018);
    addParameter(p, 'TextInterpreter', 'none');
    addParameter(p, 'FontSize', 10);
    addParameter(p, 'FontWeight', 'normal');
    parse(p, varargin{:});

    barCol = p.Results.BarColor;
    lw     = p.Results.LineWidth;
    xEarly = p.Results.EarlyX;
    xLate  = p.Results.LateX;
    toffs  = p.Results.TextOffset;
    lspace = p.Results.LineSpacing;
    interp = p.Results.TextInterpreter;
    fsz    = p.Results.FontSize;
    fwt    = p.Results.FontWeight;

    % Geometry
    m1 = mean(xEarly);
    m2 = mean(xLate);
    xCenter = mean([m1 m2]);

    ax = gca;
    wasHold = ishold(ax);
    hold(ax, 'on');

    % % Draw bars
    % h.lowerEarly = line(ax, xEarly, [group_y group_y], 'Color', barCol, 'LineWidth', lw);
    % h.lowerLate  = line(ax, xLate,  [group_y group_y], 'Color', barCol, 'LineWidth', lw);
    % h.v1 = line(ax, [m1 m1], [group_y y_top], 'Color', barCol, 'LineWidth', lw);
    % h.v2 = line(ax, [m2 m2], [group_y y_top], 'Color', barCol, 'LineWidth', lw);
    % h.top = line(ax, [m1 m2], [y_top y_top], 'Color', barCol, 'LineWidth', lw);
    
    h.lowerEarly = line(ax, xEarly, [group_y group_y], ...
    'Color', barCol, 'LineWidth', lw, ...
    'LineStyle','-', 'Marker','none');

    h.lowerLate  = line(ax, xLate, [group_y group_y], ...
        'Color', barCol, 'LineWidth', lw, ...
        'LineStyle','-', 'Marker','none');

    h.v1 = line(ax, [m1 m1], [group_y y_top], ...
        'Color', barCol, 'LineWidth', lw, ...
        'LineStyle','-', 'Marker','none');

    h.v2 = line(ax, [m2 m2], [group_y y_top], ...
        'Color', barCol, 'LineWidth', lw, ...
        'LineStyle','-', 'Marker','none');

    h.top = line(ax, [m1 m2], [y_top y_top], ...
        'Color', barCol, 'LineWidth', lw, ...
        'LineStyle','-', 'Marker','none');


    % Stack colored labels
    y0 = y_top + toffs;

    nLines = numel(sig_texts);
    assert(size(colors,1) >= nLines && size(colors,2)==3, ...
        'colors must be N×3 (N ≥ number of text lines).');

    h.txt = gobjects(nLines,1);
    for i = 1:nLines
        thisY = y0 + (i-1) * lspace;
        h.txt(i) = text(ax, xCenter, thisY, sig_texts{i}, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'Interpreter', interp, ...
            'FontSize', fsz, ...
            'FontWeight', fwt, ...
            'Color', colors(i,:), ...
            'Clipping','off');
    end

    if ~wasHold, hold(ax, 'off'); end
end
