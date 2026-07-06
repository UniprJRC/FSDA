function cmap = sequential_hcl_matlab(n, varargin)
%SEQUENTIAL_HCL_MATLAB Colorblind-safe sequential colormap via HCL space.
%   Single-hue ramp, monotonic luminance dark->light, same construction
%   family as diverging_hcl_matlab / colorspace::sequential_hcl.
p = inputParser;
addParameter(p,'H',260);     % blue, matches viridis-ish cool end; change to taste
addParameter(p,'C1',60);
addParameter(p,'L1',20);
addParameter(p,'L2',95);
addParameter(p,'Power',1.2);
parse(p,varargin{:});
o = p.Results;

t = linspace(0,1,n)';                 % 0 = dark/extreme, 1 = light
L = o.L1 + (o.L2-o.L1).*t.^o.Power;
C = o.C1.*(1-t).^o.Power;

cmap = local_hcl2rgb(o.H*ones(n,1), C, L);
cmap = max(0, min(1, cmap));
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
