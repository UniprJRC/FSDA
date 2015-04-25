function Ytra=normYJ(Y,ColtoTra,la,Jacobian)
%normYJ computes (normalized) Yeo-Johnson transformation
%
%<a href="matlab: docsearchFS('normYJ')">Link to the help function</a>
%
%  Required input arguments:
%
%         Y :   n x v data matrix; n observations
%               and v variables
%               Rows of Y represent observations, and columns represent
%               variables.
%   ColToTra:   k x 1 integer vector specifying the variables which must be
%               transformed. If it is missing and length(la)=v all
%               variables are transformed
%        la :   k x 1 vector containing set of transformation
%               parameters for the k ColtoTra.
%
% Optional input arguments:
%
%  Jacobian :   Boolean. If true (default) the transformation is normalized
%               to have Jacobian equal to 1
%
% Output:
%
%   Ytra    : n x v data matrix containing transformed observations
%             The Yeo-Johnson transformation is the Box-Cox transformation
%             of y+1 for nonnegative values, and of |y|+1 with parameter
%             2-lambda for y negative.
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
% Yeo, I.-K. and Johnson, R. (2000) A new family of power
% transformations to improve normality or symmetry. Biometrika, 87,
% 954-959.
%
%
% See also normBoxCox
%
%
%<a href="matlab: docsearchFS('normYJ')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % Transform value -3, -2, ..., 3
    y=(-3:3)';
    lambda=0
    y1=normYJ(y,1,lambda)
    plot(y,y1)
    xlabel('Original values')
    ylabel('Transformed values')

%}

%{
    % Comparison between Box-Cox and Yeo-Johnson transformation
    close all
    y=(-2:0.1:2)';
    n=length(y);
    la=-1:1:3;
    nla=length(la);
    YtraYJ=zeros(n,nla);
    YtraBC=nan(n,nla);
    posy=y>0;
    for j=1:nla
      YtraYJ(:,j)=normYJ(y,1,la(j),false);

      YtraBC(posy,j)=normBoxCox(y(posy),1,la(j),false);
    end
    subplot(1,2,1)
    plot(y,YtraYJ)
    for j=1:nla
        text(y(1), YtraYJ(1,j),['\lambda=' num2str(la(j))])
    end

    xlabel('Original values')
    ylabel('Transformed values')
    title('Yeo-Johnson transformation')

    subplot(1,2,2)
    plot(y,YtraBC)
    xlim([y(1) y(end)])
    for j=1:nla
        text(y(16), YtraBC(22,j),['\lambda=' num2str(la(j))])
    end
    xlabel('Original values')
    ylabel('Transformed values')
    title('Box-Cox transformation')
%}

%{
    % Mussels data.
    load('mussels.mat');
    Y=mussels.data;
    la=[0.5 0 0.5 0 0];
    % Transform all columns of matrix Y according to the values of la
    Y=normYJ(Y,[],la);
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


if nargin<4
    Jacobian=true;
end

%% Normalized Yeo-Johnson transformation of columns ColtoTra using la
Ytra=Y;
for j=1:length(ColtoTra);
    cj=ColtoTra(j);
    laj=la(j);
    Ycj=Y(:,cj);
    
    nonnegs = Ycj >= 0;
    negs = ~nonnegs;
    % YJ transformation is the Box-Cox transformation of
    % y+1 for nonnegative values of y
    if laj ~=0
        Ytra(nonnegs,cj)= ((Y(nonnegs,cj)+1).^laj-1)/laj;
    else
        Ytra(nonnegs,cj)= log(Y(nonnegs,cj)+1);
    end
    
    % YJ transformation is the Box-Cox transformation of
    %  |y|+1 with parameter 2-lambda for y negative.
    if 2-laj~=0
        Ytra(negs,cj) = - ((-Y(negs,cj)+1).^(2-laj)-1)/(2-laj);
    else
        Ytra(negs,cj) = -log(-Y(negs,cj)+1);
    end
    
    % If Jacobian ==true the transformation is normalized so that its
    % Jacobian will be 1
    
    if Jacobian ==true
        Ytra(:,cj)=Ytra(:,cj) * (exp(mean(log(   (1 + abs(Y(:,cj))).^(2 * nonnegs - 1)) )))^(1 - laj);
    end
end

end