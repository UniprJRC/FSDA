function [Xback] = pivotCoordInv(Xilr,varargin)
%pivotCoordInv transforms back the output of pivotCoord
%
%<a href="matlab: docsearchFS('pivotCoordInv')">Link to the help page for this function</a>
%
%  The pivot coordinates represent one-to-one mapping and this function
%  restores the original parts
%
%  Required input arguments:
%
%   Xilr :       Input data. Matrix or table.
%                Matrix or table containing of size n-by-(D-1) contained
%                units transformed into ilr coordinates
%
%
%  Optional input arguments:
%
%   norm :  the normalizing constant. Boolean.
%                If norm is true (default) sqrt((D-i)/(D-i+1))
%                is used else norm is set equal to 1.
%                Example - 'norm',false
%                Data Types - logical
%
%
%
%  Output:
%
%   Xback :      The original data restored. Array or table.
%             X is a matrix (or table depending on the input) which is
%             of size n-by-D. Note that the rows are rescaled in such a
%             way that the sum of each row is 1.
%
%
%  See also pivotCoord
%
%
% References:
%
% Filzmoser, P., Hron, K. and Templ, M. (2018), "Applied Compositional Data
% Analysis". Springer, Cham.
%
% Atkinson, A.C., Riani,M., Corbellini,A., Perrotta D., and Todorov,V.
% (2024), Applied Robust Statistics through the Monitoring Approach,
% Heidelberg: Springer Nature.
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('pivotCoordInv')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:
%
%{
    % call to CoordInv with all default options.
    X=rand(100,5);
    X=X./sum(X,2);
    Xilr=pivotCoord(X);
    Xback=pivotCoordInv(Xilr);
    assert(max(abs(X-Xback),[],'all')<1e-12)
%}

%{
    % call to CoordInv with option norm.
    X=rand(100,5);
    X=X./sum(X,2);
    nor=false;
    Xilr=pivotCoord(X,'norm',nor);
    Xback=pivotCoordInv(Xilr,'norm',nor);
%}


%% Beginning of code

% options

norm=true;

options     = struct('norm',norm);

if nargin>1
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    norm=options.norm;
end


X=-Xilr;
[n,D1]=size(X);
D =D1 + 1;

Y=zeros(n,D);
if norm == true
    Y(:,1) = -sqrt((D - 1)/D) * X(:, 1);
else
    Y(:,1)  =  X(:,1);
end

for i=2:D
    for j = 1:(i - 1)
        if norm == true
            fac=sqrt((D - j + 1) * (D - j));
        else
            fac=1;
        end

        Y(:,i) = Y(:,i) + X(:, j)/fac;
    end
end

for i=2:(D-1)
    if norm == true

        fac=sqrt((D - i)/(D - i + 1));
    else
        fac=1;
    end
    Y(:,i) = Y(:,i) - X(:,i) * fac;

end

yexp = exp(Y);
Xback = yexp./sum(yexp,2);


end

%FScategory:UTISTAT
