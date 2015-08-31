clear;

days = 0:5:35;
conc = [515,420,370,250,135,120,60,20];
temp = [29,23,27,25,20,23,23,27];

figure
bar(days,temp)
xlabel('Day', 'fontsize', 16)
ylabel('Temperature (^{o}C)', 'fontsize', 16)

% Save the handle to the axes, hAxes, and get its position.
hAxes = gca;
set(hAxes,'fontsize', 16);
hAxes_pos = get(hAxes,'Position');

% Create a second axes at the same position as the first axes. Store the handle to the second axes, hAxes2. Plot the concentration data using a thick line.
hAxes2 = axes('Position',hAxes_pos);
plot(days,conc,'LineWidth',3,'Color',[0,.7,0.7]);

% Ensure that the second axes does not interfere with the first axes by changing the y-axis location to the right side of the graph and making the background clear. 
% Remove the x-axis tick mark labels of the second axes by setting the XTickLabel property to an empty array.
set(hAxes2,'YAxisLocation','right','Color','none','XTickLabel',[], 'fontsize', 16);

% To align the x-axis of both axes, set the XLim property of the second axes equal to the XLim property of the first axes.
h1_xlim = get(hAxes,'XLim'); % store x-axis limits of first axes
set(hAxes2,'XLim',h1_xlim, 'fontsize', 16) % specify x-axis limits of second axes

% Add a title, an axis label for the right y-axis, and a text annotation.
ylabel('Concentration')
title('Trend Chart for Concentration')
text(12,380,'Concentration','Rotation',-54,'FontSize',16,'Color',[0,.7,0.7])

