function y = rescale(x,a,b)
% rescale data in [a,b]
%
%   y = rescale(x,a,b);
%
% Copyright 2008-2014.
% Written by FSDA team
%
%
%<a href="matlab: docsearch('rescale')">Link to the help function</a>
% Last modified 08-Dec-2013

%% Beginning of code
if nargin<2
    a = 0;
end
if nargin<3
    b = 1;
end

m = min(x(:));
M = max(x(:));

if M>m
    y = (b-a) * (x-m)/(M-m) + a;
else
    y = x;
end
