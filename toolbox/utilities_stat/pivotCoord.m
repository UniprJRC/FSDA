function [Xilr] = pivotCoord(X,varargin)
%pivotCoord transforms into Pivot coordinates as a special case of isometric logratio coordinates
%
%<a href="matlab: docsearchFS('pivotCoord')">Link to the help page for this function</a>
%
%  Pivot coordinates map D-part compositional data from the simplex into a
%  (D-1)-dimensional real space isometrically. From our choice of pivot
%  coordinates, all the relative information about one of parts (or about
%  two parts) is aggregated in the first coordinate (or in the first two
%  coordinates in case of symmetric pivot coordinates, respectively).
%
%  Required input arguments:
%
%   X :          Input data. Matrix or table.
%                Matrix or table containing non negative values
%
%
%  Optional input arguments:
%
%   norm :  the normalizing constant. Boolean.
%                If norm is true (default) sqrt((D-i)/(D-i+1))
%                is used else norm is set equal to 1.
%                This option takes effect just if sympivotcoord is
%                false 
%                Example - 'norm',false
%                Data Types - logical
%
%   pivotvar :  pivotal variable. Positive integer between 1 and D.
%               If any other number than 1, the data are resorted in that
%               sense that the pivotvar is shifted to the first
%               part.
%                Example - 'pivotvar',1
%                Data Types - scalar
%
% sympivotcoord :  symmetric pivot coordinate. Boolean.
%                As default symmetric pivot coordinates are not used.
%                If sympivotcoord is true the output matrix is of size n-by-2 else
%                it is of size n-by-(D-1).
%                Example - 'sympivotcoord',true
%                Data Types - logical
%
%
%  Output:
%
%   Xilr :      The data represented in pivot coordinates. Array of table.
%               Xilr is a matrix (or table depending on the input) which is
%               of size n-by-2, if sympivotcoord is true or  n-by-(D-1) if
%               sympivotcoord is false.
%
%
%  See also pivotCoordInv
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
%<a href="matlab: docsearchFS('pivotCoord')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:
%
%{
    % Call to pivotCoord with all the default options.
    % First variable as pivot variable
    X=rand(100,5);
    Xilr=pivotCoord(X);
%}

%{
    % Call to pivotCoord with options.
    % Third variable as pivot variable
    X=rand(100,5);
    Xilr=pivotCoord(X,'pivotvar',3);
%}


%{
    % Call to pivotCoord with option norm.
    % Third variable as pivot variable and no nnormalizing constant
    X=rand(100,5);
    Xilr=pivotCoord(X,'pivotvar',3,'norm',false);
%}

%{
    % Call to pivotCoord with option sympivotcoord.
    X=rand(100,5);
    % In this case the output is 100-by-2
    Xilr=pivotCoord(X,'sympivotcoord',true);
%}

%{
    %% Call to pivotCoord with input a table.
    load car
    Xt=car(:,2:6);
    % In this case the rowNames are those of the original input table
    % The variable names are z1, ..., z4
    Xilr=pivotCoord(Xt);
    head(Xilr)
%}


%% Beginning of code
sympivotcoord=false;
pivotvar=1; % the pivotal variable is assumed to be the first
norm=true;

[n,D]=size(X);

options     = struct('sympivotcoord',false,'pivotvar',pivotvar,'norm',norm);

if nargin>1
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    sympivotcoord=options.sympivotcoord;
    pivotvar=options.pivotvar;
    norm=options.norm;
end

if istable(X)
    istab=true;
    if sympivotcoord==true
        Z=X(:,1:2);
        Z.Properties.VariableNames="z"+(1:2);
    else
        Z=X(:,1:(D-1));
        Z.Properties.VariableNames="z"+(1:(D-1));
    end
    X=X{:,:};
else
    istab=false;
end


if any(X(:)<0)
    error('FSDA:pivotCoord:WrongInputOpt','negative values not allowed')
end

if sympivotcoord==true
    Xilr=zeros(n,2);
    p1=sqrt(D - 1 + sqrt(D * (D - 2)))/sqrt(2 * D);
    p2=prod(X(:,3:D), 2);
    p3 = (sqrt(D - 2) + sqrt(D))/(sqrt(D - 2) * (D - 1 + sqrt(D * (D - 2))));
    p4 = 1/(D - 1 + sqrt(D * (D - 2)));
    Xilr(:,1) = p1 * (log(X(:, 1)./(X(:,2).^p4 .* p2.^p3)));
    Xilr(:,2) = p1 * (log(X(:,2)./(X(:,1).^p4 .* p2.^p3)));
else
    X=X(:,[pivotvar setdiff(1:D,pivotvar)]);
    Xilr=zeros(n,D-1);

    for i=1:D-1
        if norm==true
            fac=sqrt((D-i)/(D-i+1));
        else
            fac=1;
        end
        Xilr(:,i)=fac*log(geomean(X(:,i+1:end),2)./X(:,i));
    end
    Xilr=-Xilr;
end

% If the input is a table the output is also a table with variable names
%z1, z2,..., zD
if istab==true
    Z{:,:}=Xilr;
    Xilr=Z;
end

%FScategory:UTISTAT
