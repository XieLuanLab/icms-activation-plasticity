function fig2svg(name, w, h, folder)
    if nargin < 4
        % Determine figure folder from name (e.g., "Fig2C" -> output/fig2/)
        tokens = regexp(name, '(?i)(S?Fig)[_]?(\d+)', 'tokens', 'once');
        if ~isempty(tokens)
            prefix = lower(tokens{1});
            num = tokens{2};
            if startsWith(prefix, 'sfig')
                subfolder = sprintf('fig_s%s', num);
            else
                subfolder = sprintf('fig%s', num);
            end
            folder = fullfile(pwd, 'output', subfolder);
        else
            folder = fullfile(pwd, 'output');
        end
    end
    if ~exist(folder,'dir'), mkdir(folder); end
    f = gcf;
    set(f, 'Units','inches', 'Position',[1 1 w h], ...
           'PaperUnits','inches', 'PaperSize',[w h], ...
           'PaperPosition',[0 0 w h], ...
           'Renderer','painters');
    file = fullfile(folder, name + ".svg");
    print(f, file, '-dsvg', '-vector');
    fprintf('  Saved: %s\n', file);
end
