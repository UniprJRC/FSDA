function Ytra=basicPower(Y,ColtoTra,la, varargin)
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
%  inverse :    Inverse transformation. Logical. If inverse is true, the
%               inverse transformation is returned. The default value of
%               inverse is false.
%                 Example - 'inverse',true
%                 Data Types - Logical
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
% Box, G.E.P. and Cox, D.R. (1964), An analysis of transformations (with
% Discussion), "Journal of the Royal Statistical Society Series B", 
% Vol. 26, pp. 211-252.
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('basicPower')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


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

%{
    % Simulated data check inverse transformation.
    n=100;p=5;
    Y=abs(randn(n,p));
    la=[0.5 0 -0.5 2 0];
    % Transform all columns of matrix Y according to the values of la
    Ytra=basicPower(Y,[],la);
    Ychk=basicPower(Ytra,[],la,'inverse',true);
    disp(max(max(abs(Y-Ychk))))
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

inverse=false;

if nargin>3
    options=struct('inverse',inverse);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:normBoxCox:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    inverse=options.inverse;
end



%% basic power transformation of columns ColtoTra using la
Ytra=Y;

if inverse== false
    
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
else
    
    for j=1:length(ColtoTra)
        cj=ColtoTra(j);
        laj=la(j);
        if  laj~=0
            Ytra(:,cj)=Y(:,cj).^(1/laj);
        else
            Ytra(:,cj)=exp(Y(:,cj));
        end
    end
end
end

%FScategory:UTISTAT