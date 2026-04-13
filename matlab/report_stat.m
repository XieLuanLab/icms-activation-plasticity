function report_stat(csvPath, panel, comparison, vals1, vals2)
% REPORT_STAT  Append one row of Mann-Whitney U stats to a CSV file.
%   report_stat(csvPath, panel, comparison, vals1, vals2)
%
%   Columns: panel, comparison, n1, n2, median1, median2, test, U, p, r

    n1 = numel(vals1);
    n2 = numel(vals2);
    med1 = median(vals1, 'omitnan');
    med2 = median(vals2, 'omitnan');

    [p, ~, stats] = ranksum(vals1, vals2);

    % Compute U from rank sum
    W = stats.ranksum;
    U = W - n1*(n1+1)/2;

    % Effect size r = |Z| / sqrt(N)
    if isfield(stats, 'zval')
        Z = stats.zval;
        r = abs(Z) / sqrt(n1 + n2);
    else
        Z = NaN;
        r = NaN;
    end

    % Write header if file doesn't exist
    if ~exist(csvPath, 'file')
        fid = fopen(csvPath, 'w');
        fprintf(fid, 'panel,comparison,n1,n2,median1,median2,test,U,p,r\n');
    else
        fid = fopen(csvPath, 'a');
    end

    fprintf(fid, '%s,%s,%d,%d,%.4f,%.4f,Mann-Whitney U,%.1f,%.4g,%.4f\n', ...
        panel, comparison, n1, n2, med1, med2, U, p, r);
    fclose(fid);

    fprintf('  %s | %s | n=(%d,%d) | median=(%.2f,%.2f) | U=%.1f | p=%.4g | r=%.4f\n', ...
        panel, comparison, n1, n2, med1, med2, U, p, r);
end
