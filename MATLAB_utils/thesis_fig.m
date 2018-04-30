function thesis_fig(handle,filename)
% Thesis figure configurator
%   configures figure identified by handle h for insertion in thesis
%   TURNS THE GRID ON
%   adjusts axes and defines page size to appropriately size whitespace
%   saves pdf with filename

grid on
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
fig = gcf;
fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];

saveas(handle,filename,'pdf')
end
