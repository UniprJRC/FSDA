function y = rescale(x,a,b)
% rescale data in [a,b]
%
%   y = rescale(x,a,b);
%
% Copyright 2008-2013.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti 
%            and Vytis Kopustinskas (2009-2010)
%
%<a href="matlab: docsearch('rescale')">Link to the help function</a>
% Last modified 02-May-2013

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
