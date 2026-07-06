function gb = plotGeobubbleDiverging(lat, lon, sizeVar, colorVar, titleStr, sizeLegend, colorLegend, varargin)

%plotGeobubbleDiverging Geobubble plot with PC1 mapped to marker size and
%   PC2 (binned into 3 classes) mapped to a diverging colorblind-safe
%   palette. This is replacing the previously duplicated logic in
%   biplotAPP.m and computePCA.m.
%
%   lat, lon      : coordinates
%   sizeVar       : continuous variable mapped to bubble size (e.g. PC1)
%   colorVar      : continuous variable to be discretized and mapped to
%                   color (e.g. PC2)
%   titleStr      : figure/gb title
%   sizeLegend    : SizeLegendTitle string
%   colorLegend   : ColorLegendTitle string
%   varargin      : extra name-value pairs forwarded to geobubble
%                   (e.g. 'MapLayout','maximized')

edges = linspace(min(colorVar), max(colorVar), 4);
cate = discretize(colorVar, edges, 'categorical');
cate = reordercats(cate, string(categories(cate)));  % ensure numeric bin order

gb = geobubble(lat, lon, sizeVar, cate, ...
    'Basemap','topographic', ...
    varargin{:});
gb.BubbleColorList   = diverging_hcl_matlab(3);
gb.Title             = titleStr;
gb.SizeLegendTitle   = sizeLegend;
gb.ColorLegendTitle  = colorLegend;
end
