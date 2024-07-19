% Example data
x = linspace(0, 10, 100);
y1 = sin(x);
y2 = 5 * cos(x);

% Create the first plot (left y-axis)
figure; % Create a new figure window
plot(x, y1, 'b-');
xlabel('X-axis');
ylabel('Y1-axis', 'Color', 'b');
title('Double Axes Plot');

% Hold the plot to overlay the second plot
hold on;

% Add second plot with a different scale (right y-axis)
plot(x, y2, 'r-');

% Adjust the axes properties and labels
ax1 = gca; % Get current axes (left y-axis)
ax1.XColor = 'k'; % Set x-axis color for the first plot

% Create a second y-axis on the right side
ax2 = axes('Position', ax1.Position, 'YAxisLocation', 'right', 'Color', 'none');
ax2.YColor = 'r'; % Set y-axis color for the second plot
ylabel('Y2-axis', 'Color', 'r');

% Add legend
legend('Y1-axis', 'Y2-axis');

% Adjust plot layout if necessary
grid on;

