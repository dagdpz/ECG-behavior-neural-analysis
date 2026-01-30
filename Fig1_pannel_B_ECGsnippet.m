% Step 1: Open figure (invisible so it doesn't pop up)
fig = openfig('Y:\Projects\_Shamim\Bacchus_bsa_ecg_results\20210720_block_04.fig');
axesHandles = findall(fig, 'Type', 'axes');
[~, idx] = sort(arrayfun(@(ax) ax.Position(2), axesHandles), 'descend'); 
axSorted = axesHandles(idx);
ax1 = axSorted(1);   % first subplot
ax2 = axSorted(2);   % second subplot% Step 4: Change x limits
ax3 = axSorted(3);

xlim(ax1, [2370 2375]);
legend(ax1, 'show', ...
       'Location', 'southoutside', ...
       'Orientation', 'horizontal');
xlim(ax2, [2370 2375]);
% xlim(ax3, [2370 2375]);
figname = fullfile(['Y:\Projects\Pulv_bodysignal\LFP\Figures',filesep,'Fig1_PanelB_ECGsnippet']);
export_fig(fig,[figname,'.pdf']);
