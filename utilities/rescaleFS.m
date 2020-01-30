function y = rescaleFS(x,a,b)
%rescale rescales numeric array to have specified minimum (a)  and maximum (b)
%
%<a href="matlab: docsearchFS('rescaleFS')">Link to the help function</a>
%
%  Required input arguments:
%     x : vector or matrix or 3D array of elements which must be rescaled
%
%  Optional input arguments:
%     a : scalar, minimum of the required interval 
%     (default value of a is 0)
%
%     b : scalar, maximum of the required interval
%     (default value of b is 1)
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('rescaleFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Rescale random numbers in the interval [3 4]
    x = abs(randn(100,1));
    % Rescale the elements of vector x  in such a way their minimum is 3
    % and their maximum is 4
    x=rescaleFS(x,3,4);
    disp(['min(x)=' num2str(min(x))])
    disp(['max(x)=' num2str(max(x))])
 
%}

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
end
