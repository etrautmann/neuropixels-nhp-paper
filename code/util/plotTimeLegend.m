function figh = plotTimeLegend(fileLength, numBins, cmap)
    % Ensure the colormap has enough colors for the specified number of bins
    cmap = colormap(cmap);  % Load the colormap
    if size(cmap, 1) < numBins
        error('The colormap does not have enough colors for the specified number of bins.');
    end
    
    % Generate bin edges and corresponding x values for each segment
    binEdges = linspace(0, fileLength, numBins + 1);
    
    % Plot each segment with its respective color
    figh = figure;
    hold on;
    for i = 1:numBins
        % Define the start and end points of the line segment
        x = [binEdges(i), binEdges(i+1)];
        y = [0, 0]; % constant y value to create a horizon line
        
        % Plot the segment with the appropriate color
        plot(x, y, 'Color', cmap(i, :), 'LineWidth', 10); % adjust LineWidth for visibility
    end
    hold off;
    
    % Formatting the plot
    xlim([0, fileLength]);
    ylim([-0.5, 0.5]);  % narrow y-axis range to keep line centered
    
    % Remove y-axis, only keep x-axis as requested
    set(gca, 'ytick', []);
    ylabel([]);
    xlabel('Time (seconds)');
    title('Horizon Line Legend');
    % axis off;
    axis tight;
    makePrettyAxis('xOnly', 'true')

    
end