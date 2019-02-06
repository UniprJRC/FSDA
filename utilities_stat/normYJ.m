function Ytra=normYJ(Y,ColtoTra,la, varargin)
%normYJ computes (normalized) Yeo-Johnson transformation
%
%<a href="matlab: docsearchFS('normYJ')">Link to the help function</a>
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
%
% Optional input arguments:
%
%  Jacobian :   Requested Jacobian of transformed values. true (default) or
%               false. If true (default) the transformation is normalized
%               to have Jacobian equal to 1
%                 Example - 'Jacobian',true
%                 Data Types - Logical
%
%  inverse :    Inverse transformation. Logical. If inverse is true, the
%               inverse transformation is returned. The default value of
%               inverse is false.
%                 Example - 'inverse',true
%                 Data Types - Logical
%
% Output:
%
%   Ytra    : transformed data matrix. Matrix. n x v data matrix containing
%               transformed observations
%             The Yeo-Johnson transformation is the Box-Cox transformation
%             of y+1 for nonnegative values, and of |y|+1 with parameter
%             2-lambda for y negative.
%
%
% See also normBoxCox
%
% References:
%
% Yeo, I.K and Johnson, R. (2000), A new family of power transformations to
% improve normality or symmetry, "Biometrika", Vol. 87, pp. 954-959.
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('normYJ')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Example of use of normYJ with all default options.
    % Transform value -3, -2, ..., 3
    y=(-3:3)';
    lambda=0
    y1=normYJ(y,1,lambda)
    plot(y,y1)
    xlabel('Original values')
    ylabel('Transformed values')

%}

%{
    % Comparison between Box-Cox and Yeo-Johnson transformation.
    close all
    y=(-2:0.1:2)';
    n=length(y);
    la=-1:1:3;
    nla=length(la);
    YtraYJ=zeros(n,nla);
    YtraBC=nan(n,nla);
    posy=y>0;
    for j=1:nla
      YtraYJ(:,j)=normYJ(y,1,la(j),'Jacobian',false);

      YtraBC(posy,j)=normBoxCox(y(posy),1,la(j),'Jacobian',false);
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
    % Simulated data check inverse transformation.
    n=100;p=5;
    Y=randn(n,p);
    Y(3,1:3)=0;
    la=[0.5 0 -0.5 2 0];
    % Transform all columns of matrix Y according to the values of la
    Ytra=normYJ(Y,[],la,'Jacobian',false);
    Ychk=normYJ(Ytra,[],la,'Jacobian',false,'inverse',true);
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

Jacobian=true;
inverse=false;

if nargin>2
    options=struct('Jacobian',Jacobian,'inverse',inverse);
    
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
    Jacobian=options.Jacobian;
    inverse=options.inverse;
end

%% Normalized Yeo-Johnson transformation of columns ColtoTra using la
Ytra=Y;

if inverse== false
    
    for j=1:length(ColtoTra)
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
    
else % inverse transformation
     for j=1:length(ColtoTra)
        cj=ColtoTra(j);
        laj=la(j);
        Ycj=Y(:,cj);
        
        nonnegs = Ycj >= 0;
        negs = ~nonnegs;
        
        
        % YJ transformation is the Box-Cox transformation of
        % y+1 for nonnegative values of y
        if laj ~=0
            Ytra(nonnegs,cj)=  (laj*Y(nonnegs,cj)+1).^(1/laj)   -1 ;
        else
            Ytra(nonnegs,cj)= expm1(Y(nonnegs,cj));
        end
        
        if 2-laj~=0
            Ytra(negs,cj) = 1 - (    ((laj-2)*Y(negs,cj)+1).^(1/(2-laj))  ) ;
        else
            Ytra(negs,cj) = -expm1(-Y(negs,cj));
        end
        
     end
     
    % insert a NaN every time there is a number which is not real
    Ytra(imag(Ytra(:))~=0)=NaN;

%             ans[index] <- (y[index] * lambda[index] + 1)^(1/lambda[index]) - 
%                 1
%         if (any(index <- y >= 0 & abs(lambda) <= epsilon)) 
%             ans[index] <- expm1(y[index])
%         if (any(index <- y < 0 & abs(lambda - 2) > epsilon)) 
%             ans[index] <- 1 - (-(2 - lambda[index]) * y[index] + 
%                 1)^(1/(2 - lambda[index]))
%         if (any(index <- y < 0 & abs(lambda - 2) <= epsilon)) 
%             ans[index] <- -expm1(-y[index])  
end

end
%FScategory:UTISTAT