function [negColor,posColor] = plotLoadingsDiverging(loadings, varnames, NumComponents, figName, figTag)
%plotLoadingsDiverging Plots PCA loadings bar charts with a diverging
%   colorblind-safe palette (sign-encoded), replacing the previous
%   duplicated flat-green bar plot logic that lived separately in
%   biplotAPP.m and computePCA.m.
%
%   loadings       : n_vars x NumComponents matrix of loadings
%   varnames       : cellstr/string array of variable names
%   NumComponents  : number of PCs to plot (one subplot each)
%   figName        : figure 'Name' property (e.g. 'Loadings')
%   figTag         : figure 'Tag' property (e.g. 'pl_loadings')

delete(findobj(0, 'type', 'figure','tag',figTag));
figure('Name',figName,'Tag',figTag)

xlabels = categorical(varnames,varnames);

cmapDiv  = aux.diverging_hcl_matlab(256);
negColor = cmapDiv(1,:);
posColor = cmapDiv(end,:);

for i = 1:NumComponents
    subplot(NumComponents,1,i)
    vals = loadings(:,i);

    cdata = repmat(posColor, numel(vals), 1);
    cdata(vals < 0, :) = repmat(negColor, sum(vals<0), 1);

    b = bar(xlabels, vals);
    b.FaceColor = 'flat';
    b.CData = cdata;

    xtips = b(1).XData;
    ytips = b(1).YData;
    barlabels = string(round(vals,2));
    text(xtips,ytips,barlabels,'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom')
    title(['Correlations  with PC' num2str(i)])
end
end
