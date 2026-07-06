function cmap = diverging_hcl_matlab(n, varargin)

%DIVERGING_HCL_MATLAB Colorblind-safe diverging colormap via HCL space.
%
% Blue-Red style diverging HCL palette, following Zeileis et al. (2020).
% This palette is used for a signed quantity with a meaningful
% zero, such as PCA loadings or correlation coefficients, for which we
% need a colorblind-safe diverging red-blue palette centered at zero,
% rather than a sequential or single-hue scale that would obscure the
% sign of the values.
%
% Example of use:
% cmapBackground = diverging_hcl_matlab(256); 
% colormap(cmapBackground);
%
%   cmap = DIVERGING_HCL_MATLAB(n) returns an n-by-3 RGB colormap that
%   interpolates from a blue extreme, through a neutral (near-white)
%   center, to a red extreme, following the same HCL (hue-chroma-
%   luminance) construction used by the R package colorspace
%   (Zeileis et al., 2020, JSS) for diverging_hcl().
%
%   Unlike RGB-linear interpolation (e.g. colorRampPalette in base R),
%   this varies luminance and chroma monotonically and symmetrically in
%   a perceptually uniform space, so that equal steps in the data (and
%   in the colormap index) correspond to approximately equal perceived
%   color differences on both sides of the neutral midpoint.
%
%   cmap = DIVERGING_HCL_MATLAB(n, 'Name', Value, ...) allows overriding
%   defaults:
%       'H1'    hue of the low (negative) extreme   [default 260, blue]
%       'H2'    hue of the high (positive) extreme  [default 12,  red]
%       'C1'    max chroma at the extremes          [default 80]
%       'L1'    luminance at the extremes           [default 30]
%       'L2'    luminance at the center              [default 95]
%       'Power' easing exponent for chroma/luminance [default 1.5]
%
%   Example:
%       cmap = diverging_hcl_matlab(256);
%       colormap(cmap); colorbar; caxis([-1 1]);

p = inputParser;
addParameter(p,'H1',260);
addParameter(p,'H2',12);
addParameter(p,'C1',80);
addParameter(p,'L1',30);
addParameter(p,'L2',95);
addParameter(p,'Power',1.5);
parse(p,varargin{:});
o = p.Results;

half = ceil(n/2);
t = linspace(1,0,half)';              % 1 at extreme, 0 at center
L = o.L2 - (o.L2-o.L1).*t.^o.Power;
C = o.C1.*t.^o.Power;

negHalf = local_hcl2rgb(o.H1*ones(half,1), C, L);   % extreme -> center
posHalf = local_hcl2rgb(o.H2*ones(half,1), C, L);   % extreme -> center

cmap = [negHalf; flipud(posHalf(2:end,:))];         % extreme -> center -> extreme
cmap = max(0, min(1, cmap));

if size(cmap,1) ~= n
    idx  = round(linspace(1, size(cmap,1), n));
    cmap = cmap(idx,:);
end
end

function rgb = local_hcl2rgb(H, C, L)

% Polar LUV (H,C,L) -> LUV -> CIEXYZ -> linear sRGB -> gamma-corrected sRGB
Hr = deg2rad(H);
U = C.*cos(Hr);
V = C.*sin(Hr);

% D65 white point in CIE 1976 u',v'
Un = 0.19783982; Vn = 0.46833630; Yn = 1;

Y = Yn.*((L+16)/116).^3;
smallL = L <= 8;
Y(smallL) = Yn.*L(smallL)/903.3;

u0 = U./(13*L) + Un;
v0 = V./(13*L) + Vn;

X = Y.*9.*u0./(4*v0);
Z = Y.*(12 - 3*u0 - 20*v0)./(4*v0);

M = [ 3.2404542 -1.5371385 -0.4985314;
     -0.9692660  1.8760108  0.0415560;
      0.0556434 -0.2040259  1.0572252];
rgb = [X Y Z]*M';
rgb = real(rgb);

mask = rgb > 0.0031308;
rgb(mask)  = 1.055*rgb(mask).^(1/2.4) - 0.055;
rgb(~mask) = 12.92*rgb(~mask);
end
