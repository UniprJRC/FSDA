function Ytra=normYJpn(Y,ColtoTra,la, varargin)
%normYJ computes (normalized) extended Yeo-Johnson transformation
%
%<a href="matlab: docsearchFS('normYJpn')">Link to the help function</a>
%
% The transformations for negative and positive responses were determined
% by Yeo and Johnson (2000) by imposing the smoothness condition that the
% second derivative of zYJ(λ) with respect to y be smooth at y = 0. However
% some authors, for example Weisberg (2005), query the physical
% interpretability of this constraint which is oftern violated in data
% analysis. Accordingly, Atkinson et al (2019) and (2020) extend the
% Yeo-Johnson transformation to allow two values of the transformations
% parameter: λN for negative observations and λP for non-negative ones.
%
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
%        la :   transformation parameters. Matrix.
%               k x 2 vector containing set of transformation
%               parameters for the k ColtoTra. The first column contains
%               the transformation parameter for positive observations and
%               the second column the transformation parameter for negative
%               observations.
%                Data Types - single|double
%
%
% Optional input arguments:
%
%  inverse :    Inverse transformation. Logical. If inverse is true, the
%               inverse transformation is returned. The default value of
%               inverse is false.
%                 Example - 'inverse',true
%                 Data Types - Logical
%
%  Jacobian :   Requested Jacobian of transformed values. true (default) or
%               false. If true (default) the transformation is normalized
%               to have Jacobian equal to 1
%                 Example - 'Jacobian',true
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
% See also normBoxCox, normYJ
%
% References:
%
% Atkinson, A.C. Riani, M., Corbellini A. (2019), The analysis of
% transformations for profit-and-loss data, Journal of the Royal
% Statistical Society, Series C, "Applied Statistics",
% https://doi.org/10.1111/rssc.12389
% Atkinson, A.C. Riani, M. and Corbellini A. (2020), The Box-Cox
% Transformation: Review and Extensions, "Statistical Science", in press.
% Yeo, I.K and Johnson, R. (2000), A new family of power transformations to
% improve normality or symmetry, "Biometrika", Vol. 87, pp. 954-959.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('normYJpn')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Example of use of normYJext with all default options.
    % Transform value -3, -2, ..., 3
    y=(-3:3)';
    lambda=[-0.5 0.5]
    y1=normYJpn(y,1,lambda)
    plot(y,y1)
    xlabel('Original values')
    ylabel('Transformed values')

%}

%{
    % Compare Yeo & Johnson with extended Yeo and Yohnson.
    % Transform value -3, -2, ..., 3
    k=3;
    y=(-k:0.01:k)'; 
    % Two values of lambda for extended Yeo and Johnson
    lambda=[0 0.5];
    Jacobian=false;
    % Just one value of lambda for traditional Yao and Johnson
    y1=normYJ(y,1,lambda(1),'Jacobian',Jacobian);
    
    ypn=normYJpn(y,1,lambda,'Jacobian',Jacobian);
    plot(y,[y1 ypn])
    xlabel('Original values')
    ylabel('Transformed values')
    legend({['YJ with $\lambda$=' num2str(lambda(1))],...
        ['$YJ_{ext}$ with $\lambda_P$=' num2str(lambda(1)) ' and  $\lambda_N$=' num2str(lambda(2))]},...
        'Interpreter','latex') 
%}

%{
    % Simulated data check inverse transformation.
    n=100;p=5;
    Y=randn(n,p);
    Y(3,1:3)=0;
    % Different values of transformation parameters for positive and
    % negative observations
    la=[[0.5; 0; -0.5; 2; 0] , [0.5; 0; -0.5; 2; 0]+0.8];
    % Transform all columns of matrix Y according to the values of la
    Ytra=normYJpn(Y,[],la,'Jacobian',false);
    Ychk=normYJpn(Ytra,[],la,'Jacobian',false,'inverse',true);
    disp(max(max(abs(Y-Ychk))))
%}


%% Beginning of code

% Input parameters checking
% Extract size of the data
v=size(Y,2);

if nargin<1
    error('FSDA:normYJpn:missingInputs','Input data matrix is missing');
end

if nargin<2
    error('FSDA:normYJpn:missingInputs','Vector ColtoTra which specifies which variables to transform is missing');
end

if nargin<3
    error('FSDA:normYJpn:missingInputs','Vector la which specifies how to transforme the variables is missing');
end

if isempty(ColtoTra) && size(la,1)==v
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
            error('FSDA:normYJpn:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
        lajPos=la(j,1);
        lajNeg=la(j,2);
        Ycj=Y(:,cj);
        
        nonnegs = Ycj >= 0;
        negs = ~nonnegs;
        
        
        % YJ transformation is the Box-Cox transformation of
        % y+1 for nonnegative values of y
        if lajPos ~=0
            Ytra(nonnegs,cj)= ((Y(nonnegs,cj)+1).^lajPos-1)/lajPos;
        else
            Ytra(nonnegs,cj)= log(Y(nonnegs,cj)+1);
        end
        
        % YJ transformation is the Box-Cox transformation of
        %  |y|+1 with parameter 2-lambda for y negative.
        if 2-lajNeg~=0
            Ytra(negs,cj) = - ((-Y(negs,cj)+1).^(2-lajNeg)-1)/(2-lajNeg);
        else
            Ytra(negs,cj) = -log(-Y(negs,cj)+1);
        end
        
        % If Jacobian ==true the transformation is normalized so that its
        % Jacobian will be 1
        if Jacobian ==true
            % logJ=log of the Jacobian (J)
            logJ=(lajPos-1)*sum(log(Y(nonnegs,cj)+1)) ...
                - (lajNeg-1)*sum( log(-Y(negs,cj)+1) );
            % J1n=(J)^(1/n);
            J1n=(exp(logJ))^(1/size(Y,1));
            Ytra(:,cj)=Ytra(:,cj) /J1n ;
            % Ytra(:,cj)=Ytra(:,cj) * (exp(mean(log(   (1 + abs(Y(:,cj))).^(2 * nonnegs - 1)) )))^(1 - laj);
        end
    end
    
else % inverse transformation
    for j=1:length(ColtoTra)
        cj=ColtoTra(j);
        lajPos=la(j,1);
        lajNeg=la(j,2);
        
        Ycj=Y(:,cj);
        
        nonnegs = Ycj >= 0;
        negs = ~nonnegs;
        
        
        % YJ transformation is the Box-Cox transformation of
        % y+1 for nonnegative values of y
        if lajPos ~=0
            Ytra(nonnegs,cj)=  (lajPos*Y(nonnegs,cj)+1).^(1/lajPos)   -1 ;
        else
            Ytra(nonnegs,cj)= expm1(Y(nonnegs,cj));
        end
        
        if 2-lajNeg~=0
            Ytra(negs,cj) = 1 - (    ((lajNeg-2)*Y(negs,cj)+1).^(1/(2-lajNeg))  ) ;
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