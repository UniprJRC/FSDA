function Ytra=basicPower(Y,ColtoTra,la)
%basicPower computes the basic power transformation
%
%<a href="matlab: docsearchFS('basicPower')">Link to the help function</a>
%
%  Required input arguments:
%
% Y :           Input data. Matrix. 
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%   ColtoTra:   Variable to transform. Vector.  k x 1 integer vector
%               specifying the variables which must be
%               transformed. If it is missing and length(la)=v all
%               variables are transformed
%                Data Types - single|double
%        la :   transformation parameters. Vector.
%               k x 1 vector containing set of transformation
%               parameters for the k ColtoTra.
%                Data Types - single|double
%
% Optional input arguments:
%
% Output:
%
%   Ytra    : transformed data matrix. Matrix. n x v data matrix containing
%               transformed observations
%
%             When $\lambda \ne 0$
%             \[
%               ytra = y^\lambda
%             \]
%             else
%             \[ 
%               ytra = log(y) 
%             \]
%
%
% See also: normBoxCox, normYJ
% 
% References:
%
%Box, G. E. P. & Cox, D. R. (1964). An analysis of transformations (with
%Discussion). J. R. Statist. Soc. B 26, 211-252
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('basicPower')">Link to the help function</a>
% Last modified 14-06-2016


% Examples:

%{
    % Example of transformation.
    % Mussels data.
    load('mussels.mat');
    Y=mussels.data;
    la=[0.5 0 0.5 0 0];
    % Transform all columns of matrix Y according to the values of la using
    % basic power transformation
    Y1=basicPower(Y,[],la);
%}



%% Input parameters checking
% Extract size of the data
v=size(Y,2);

if nargin<1
    error('FSDA:normBoxCox:missingInputs','Input data matrix is missing');
end

if nargin<2
    error('FSDA:normBoxCox:missingInputs','Vector ColtoTra which specifies which variables to transform is missing');
end

if nargin<3
    error('FSDA:normBoxCox:missingInputs','Vector la which specifies how to transforme the variables is missing');
end

if isempty(ColtoTra) && length(la)==v
    ColtoTra=1:v;
end


%% basic power transformation of columns ColtoTra using la
Ytra=Y;
for j=1:length(ColtoTra)
    cj=ColtoTra(j);
    laj=la(j);
    Ycj=Y(:,cj);
    
    if min(Ycj)<=0 && laj==1
        % if min(Ycj)<=0 and la(cj)=1 then variable is not transformed
    elseif min(Ycj)<=0 && cj ~=1
        error('FSDA:normBoxCox:Wronglaj',['lambda=' num2str(laj) ' for column ' num2str(cj) ' but min(Ycj)=' num2str(min(Ycj))])
    else
        if laj~=0
            Ytra(:,cj)=Y(:,cj).^laj;
        else
            Ytra(:,cj)=log(Y(:,cj));
        end
    end
end

end
%FScategory:UTISTAT