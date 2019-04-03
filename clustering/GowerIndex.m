function S = GowerIndex(Y, varargin)
%GowerIndex computes matrix of similarity indexes using Gower metric
%
%
%<a href="matlab: docsearchFS('GowerIndex')">Link to the help page for this function</a>
%
%   This function uses the alternating conditional expectations algorithm
%   to find the transformations of y and X that maximise the proportion of
%   variation in y explained by X. When X is a matrix, it is transformed so
%   that its columns are equally weighted when predicting y.
%
%
% Required input arguments:
%
% Y :           Input data. 2D array or MATLAB table .
%               n x p data matrix; n observations and p variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
% Optional input arguments:
%
%       l :  type of variable. Vector. Vector of length p which
%           specifies the type of variable for each column of the input
%           matrix Y.
%           $l(j)=1$ => $j$-th variable assumes orderable values (quantitative
%           variable).
%           $l(j)=2$ => $j$-th variable is binary (just values 0 and 1).
%           $l(j)=3$ => $j$-th variable assumes categorical (unorderable) values.
%           $j =1, 2, \ldots, p$.
%           If l is omitted or is empty the routine automatically tries to detect the
%           type of variable. For example, if a variable assumes just
%           values 0 and 1  the corresponding variable is classified as
%           binary. If a column of the input dataset assumes just integer
%           values (without decimal points) and the number of unique
%           elements is not greater than 20 the corresponding variable is
%           classified as categorical,
%           else the variable is classified as quantitative.
%           Example - 'l',[3 3 1]
%           Data Types - double
%
% Output:
%
%   S   :  matrix with Gower similarity coefficients. n-by-n symmetric matrix.
%          n-by-n natrix whose i-th j-th entry contains the Gower
%          similarity index between row i and row j of input matrix Y.
%
%
%
% More About:
%
% A very popular metric for mixtures of quantitative, multistate
% categorical, and binary variables is the one based on Gower's general
% similarity coefficient (Gower, 1971) that substantially unifies Jaccard's
% coefficient (binary variables), the simple matching
% coefficient (multistate categorical variables) and normalized city block
% distance (quantitative variables).
% More specifically, given two $p$-dimensional vectors $z_{i}$ and $z_{j}$, Gower's
% similarity coefficient is defined as
% \[
%   s_{ij}=\frac{\sum_{h=1}^{p_{1}}\left(1-|z_{ih}-z_{jh}|/R_{h}\right)+a+\alpha}
%   {p_{1}+(p_{2}-d)+p_{3}}, \quad 0\leq s_{ij}\leq 1,
% \]
% where $p=p_{1}+p_{2}+p_{3}$, \(p_{1}\) is the number of continuous variables, $a$ and $d$ are the number of positive and negative matches, respectively, for the \(p_{2}\) binary variables, \(\alpha\) is the number of matches for the \(p_{3}\) multi-state categorical variables, and \(R_{h}\) is the range of the \(h\)-th continuous variable.
%
%
%
% See also: pdist, tclust
%
% References:
%
% Gower, J. C. (1971), "A general coefficient of similarity and some of its
% properties", Biometrics, pp. 857-871.
% Grane', A., and Romera R. (2018), "On Visualizing Mixed-Type Data: A
% Joint Metric Approach to Profile Construction and Outlier Detection",
% Sociological Methods & Research, Vol. 47, pp. 207-239
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('GowerIndex')">Link to the help function</a>
%
%$LastChangedDate:: 2019-03-05 16:31:26 #$: Date of the last commit

% Examples:

%{
    % Example of matrix of Gower indexes using all default options.
    Y=randn(10,4);
    S=GowerIndex(Y);
%}


%{
    % Example of use of option l.
    p=3;
    n=30;
    Y=randi(20,n,p);
    % Specify that all variables are quantitative.
    l=ones(p,1);
    S=GowerIndex(Y,'l',l);
%}


%{
    % Example 2 of use of option l.
    p=3;
    n=50;
    Y=randi(120,n,p);
    % Specify that first variable is quantitative and the other 2 are categorical.
    l=[1 3 3];
    S=GowerIndex(Y,'l',l);
%}

%% Beginning of code

l=[];
if nargin>1
    options=struct('l',l);
    
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:GowerIndex:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    l   = options.l;
end

%           l(j)=1 => j-th variable assumes orderable values (quantitative
%           variable).
%           l(j)=2 => j-th variable is binary (just values 0 and 1)
%           l(j)=3 => j-th variable assumes categorical (unorderable) values.

[n,p]=size(Y);

if isempty(l)
    booquant=false(p,1); % Initialize Boolean for quantitative variables
    boobin=booquant;  % Initialize Boolean for binary variables
    boocat=booquant; % Initialize Boolean for categorical variables
    for j=1:p
        uniYj=unique(Y(:,j));
        m=length(uniYj);
        if m==1
            error('FSDA:GowerIndex:ConstCol','A column of the dataset is made up of constant values')
        end
        
        if m==2 && max(uniYj)==1 && min(uniYj)==0 % binary variable
            boobin(j)=true;
        elseif m<=20 &&  max(uniYj-round(uniYj))==0 % categorical variable
            boocat(j)=true;
        else % quantitative variable
            booquant(j)=true;
        end
    end
else
    
    % Check that the length of supplied vector l is equal to p
    if length(l)~=p
        error('FSDA:GowerIndex:WrongInput','The length of supplied input vector l must be equal to the columns of input matrix Y')
    end
    
    % Create boolean vectors containing true in position j if the corresponding
    % j-th variable matches the specification
    booquant=l==1; % Boolean for quantitative variables
    boobin=l==2;  % Boolean for binary variables
    boocat=l==3;  % Boolean for categorical variables
end

p1=sum(booquant);


% Xcont = matrix of quantitative variables
Xcont=Y(:,booquant);

% Xbin = matrix of binary variables
Xbin=Y(:,boobin);

% Xcat = matrix of categorical variables
Xcat=Y(:,boocat);

% Similarity matrix for continuous variables
rango=max(Xcont)-min(Xcont);
C=zeros(n,n);
for i=1:n
    C(i,i)=p1;
    for j=1:i-1
        C(i,j)=p1-sum(abs(Xcont(i,:)-Xcont(j,:))./rango);
        C(j,i)=C(i,j);
    end
end

% Computation of matrix A for binary variables
A=Xbin*Xbin';

J=ones(size(Xbin));
D=(J-Xbin)*(J-Xbin)';


Alpha=zeros(n,n);
for i=1:n
    for j=i:n
        Alpha(i,j)=sum(Xcat(i,:)==Xcat(j,:));
        Alpha(j,i)=Alpha(i,j);
    end
end
% Computation of Gower similarity index
S=(C+A+Alpha)./(p-D);
end
%FScategory:CLUS-RobClaMULT