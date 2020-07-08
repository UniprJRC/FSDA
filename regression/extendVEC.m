function [ext] = extendVEC(x, index, value)
% extendVEC extends vector x inserting value in the positions in index
%
%<a href="matlab: docsearchFS('extendVEC')">Link to the help function</a>
%
% Required input arguments:
%
% x:        Vector to be extended. It can be either a row or a column vector.
% index:    A vector containing the positions where to add value inside x.
%           The index values are sorted so that the code starts to insert the element "value"
%           from the smallest position.
% value:    The element to insert in x for the extension. It could be a number, string, NaN, Inf, etc...
%
%
% Optional input arguments:
%
% Output:
%
% ext:  The extension of vector x with the element value inserted in the positions in index.
%       For index > length(x), the code inserts index-length(x) elements equal to value in the tail of x.
%       It can be either a row or a column vector, depending on input x.
%
%
% Copyright 2008-2020.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('extendVEC')">Link to the help function</a>
%
%$LastChangedDate:: 2020-06-09 16:58:50 #$: Date of the last commit
%
% Example:
%{
%% Example 1: Row vector extended with the inclusion of multiple elements "Inf" in the tail.
  
  x = [1 2 3 4 5 6];
  index = [2 5 11];
  ext = extendVEC(x, index, Inf)

%}


%% Beginning of code

% The first argument which is passed is x.
if nargin<1 || isempty(x)==1
    error('FSDA:extendVEC:MissingInputs','Input vector x not specified.');
end

[n1, n2] = size(x);

if min(n1,n2)>1
    error('FSDA:extendVEC:Wrongx','x must be a vector.');
end

% The second argument which is passed is index.
if nargin<2 || isempty(index)==1
    error('FSDA:extendVEC:MissingInputs','Input vector index not specified.');
end

[k1, k2] = size(index);

if min(k1,k2)>1
    error('FSDA:extendVEC:Wrongindex','index must be a vector.');
end

% The third argument which is passed is value.
if nargin<3 || isempty(value)==1
    error('FSDA:extendVEC:MissingInputs','Input value not specified.');
end

[k1, k2] = size(value);

if k1~=1 || k2~=1
    error('FSDA:extendVEC:Wrongvalue','value must contain a singular entry.');
end

index = sort(index);

ext = x;

if n1==1 % x is a row vector.
    for r = index
        if r == 1
            ext = [value ext];
        elseif r <= length(ext)
            ext = [ext(1:(r-1)) value ext(r:end)];
        else
            ext = [ext repmat(value, 1, r-length(ext))];
        end
    end
else % x is a column vector.
    for r = index
        if r == 1
            ext = [value; ext];
        elseif r <= length(ext)
            ext = [ext(1:(r-1)); value; ext(r:end)];
        else
            ext = [ext; repmat(value, r-length(ext), 1)];
        end
    end
end

end