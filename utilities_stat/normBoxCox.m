function Ytra=normBoxCox(Y,ColtoTra,la,varargin)
%normBoxCox computes (normalized) Box-Cox transformation
%
%<a href="matlab: docsearchFS('normBoxCox')">Link to the help function</a>
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
%             When $\lambda \ne 0$
%             if jacobian=true:
%             $ytra = (y^\lambda-1)/ (G^{(\lambda-1)} \lambda)$;
%             else if jacobian=false:
%             $ytra = (y^\lambda-1)/ \lambda$;
%             where $G$ is the geometric mean of the observations.
%             When $\lambda = 0$
%             if jacobian=true:
%             $ytra = G log(y)$;
%             else if jacobian=false:
%             $ytra = log(y)$;
%             where $G$ is the geometric mean of the observations.
%
% See also normYJ
%
%
% References:
%
% Box, G.E.P. and Cox, D.R. (1964), An analysis of transformations (with
% Discussion), "Journal of the Royal Statistical Society Series B", 
% Vol. 26, pp. 211-252.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('normBoxCox')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


% Examples:

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
    % Mussels data.
    load('mussels.mat');
    Y=mussels.data;
    la=[0.5 0 0.5 0 0];
    % Transform all columns of matrix Y according to the values of la
    Y=normBoxCox(Y,[],la);
%}

%{
    % Check the inverse tranformation.
    load('mussels.mat');
    Y=mussels.data;
    la=[0.5 0 0.5 0 0];
    % Transform all columns of matrix Y according to the values of la
    Ytra=normBoxCox(Y,[],la,'Jacobian',false);
    Ychk=normBoxCox(Ytra,[],la,'inverse',true);
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


%% Normalized Box Cox transformation of columns ColtoTra using la
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
            % If Jacobian ==true the transformation is normalized so that its
            % Jacobian will be 1
            if Jacobian ==true
                Gj=exp(mean(log(Ycj)));
            else
                Gj=1;
            end
            
            if laj~=0
                Ytra(:,cj)=(Y(:,cj).^laj-1)/(laj*(Gj^(laj-1)));
            else
                Ytra(:,cj)=Gj*log(Y(:,cj));
            end
        end
    end
else
    
    for j=1:length(ColtoTra)
        cj=ColtoTra(j);
        laj=la(j);
        
        if laj~=0
            Ytra(:,cj)=(Y(:,cj)*laj+1).^(1/laj);
        else
            Ytra(:,cj)=exp(Y(:,cj));
        end
    end
end
end

%FScategory:UTISTAT