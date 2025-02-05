function [] = plot_drift_traces3(dshift, tvec, varargin)
    % PLOT_DRIFT_TRACES
    % dshift - data matrix where each row is a sample and each column is a channel/block
    % tvec - time vector
    % varargin - optional parameter-value pairs (e.g., 'YLim', [-50 50], 'Estimator', 'mean')

    % Set up input parser
    p = inputParser;
    addRequired(p, 'dshift', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'tvec', @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'YLim', [-25, 25], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'Estimator', '', @(x) ischar(x) && any(strcmp(x, {'mean', 'median', ''})));

    % Parse inputs
    parse(p, dshift, tvec, varargin{:});
    y_limits = p.Results.YLim;
    estimator = p.Results.Estimator;

    % Ensure tvec is a column vector
    if isrow(tvec)
        tvec = tvec';
    end

    % Check if mean ± SD, median with IQR, or individual traces should be plotted
    if strcmp(estimator, 'mean')
        % Calculate mean and standard deviation across columns (channels) for each sample (row)
        trace = mean(dshift, 2);  % Mean across columns for each row
        std_trace = std(dshift, 0, 2); % Standard deviation across columns for each row
        
        % Plot mean with shaded error (mean ± SD)
        fill([tvec; flipud(tvec)], [trace + std_trace; flipud(trace - std_trace)], ...
             [0.8 0.8 1], 'EdgeColor', 'none'); % Light blue fill for error
        hold on;
        plot(tvec, trace, 'b', 'LineWidth', 2); % Plot mean line in blue
        legend('Mean ± SD', 'Mean Trace');
        
    elseif strcmp(estimator, 'median')
        % Calculate median and IQR across columns (channels) for each sample (row)
        trace = median(dshift, 2);                    % Median across columns for each row
        q25 = prctile(dshift, 25, 2);                 % 25th percentile across columns for each row
        q75 = prctile(dshift, 75, 2);                 % 75th percentile across columns for each row
        
        % Plot median with shaded IQR (25th to 75th percentile)
        fill([tvec; flipud(tvec)], [q75; flipud(q25)], ...
             [0.8 0.8 1], 'EdgeColor', 'none'); % Light blue fill for IQR
        hold on;
        plot(tvec, trace, 'r', 'LineWidth', 2); % Plot median line in red
        legend('IQR (25th - 75th percentile)', 'Median Trace');
        
    else
        % Plot individual traces
        nLines = size(dshift, 2);
        phandle = plot(tvec, fliplr(dshift), 'linewidth', 1);
        grid on;
        
        cmap = (cbrewer('div', 'RdYlBu', nLines));
        legendStr = cell(1, nLines);
        
        for ii = 1:nLines
            set(phandle(ii), 'color', cmap(ii, :));
            legendStr{ii} = sprintf('Ch. block %d', ii);
        end
        
        le = legend(legendStr);
        le.FontSize = 10;
    end

    % Set y-limits, plot labels, and appearance
    ylim(y_limits);
    xlabel('Time (s)', 'FontSize', 20);
    ylabel('Shift (µm)', 'FontSize', 20);
    % title('Probe drift', 'FontSize', 20);

    ax = gca;
    ax.FontSize = 20;

    set(gcf, 'WindowStyle', 'normal');
    set(gcf, 'position', [0 0 600 400]);
    set(gcf, 'color', 'w');
end