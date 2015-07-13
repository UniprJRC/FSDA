function thresh=RobRegrSize(n,p,robest,rhofunc,bdp,eff,sizesim,Tallis)
%RobRegrSize provides proper threshold for robust estimators to obtain an empirical size equal to 1 per cent nominal size
%
%<a href="matlab: docsearchFS('robregrsize')">Link to the help function</a>
%
%
%  Required input arguments:
%
%           n : sample size. Scalar integer.
%               Number of units of the regression dataset. 
%               REMARK - simulations have been done for
%               n=50, 60, 70, 80, 90, 100, 150, 200, 250, 300, 400, 500.
%               For other values of n the threhold are found by
%               interpolation using the two closest values smaller or
%               greater than the one which has been considered
%               Data Types -single | double | int8 | int16 | int32 | int64 |uint8 | uint16 | uint32 | uint64
%           p : number of variables. Scalar integer. Number of explanatory variables.
%               REMARK - simulations have been done for p=2, 3, ..., 10. If
%               the user supplies a value of p greater than 10 the
%               correction factors are extrapolated by fitting a simple
%               quadratic model in p.
%               Data Types -single | double | int8 | int16 | int32 | int64 |uint8 | uint16 | uint32 | uint64
%       robest: robust estimator. String. String which identifies the robust estimator which is used
%               Possibile values are:
%                 'S'    S estimators;
%                 'MM'   MM estimators;
%                 'LTS'  Least trimmed squares estimator;
%                 'LTSr' Least trimmed squares estimator reweighted.
%               If robest is missing MM estimator is used
%               Data Types -char
%      rhofunc: Weight function. String.
%               String which identifies the weight function which has been
%               used for S or MM.
%               Possibile values are
%                'TB', for Tukey biweight rho function;
%                'HA', for Hampel rho function;
%                'HY', for hyperbolic rho function;
%                'OP', for optimal rho function;
%                'ST'  Soft trimming estimator (in this case an average
%                        threshold based on the TB,HY,HA and OP is used)
%                REMARK - this value is ignored if robest is LTS or LTSr
%               If rhofunc is missing and robest is 'S' or 'MM', the
%               default value of rhofunc is 'ST'.
%               Data Types -char
%         bdp :  breakdown point. Scalar.
%               Scalar between 0 and 0.5. If robest is S, LTS or LTSr
%               and bdp is missing a value of 0.5 is used as default.
%               REMARK - simulations have been done for bdp=0.25 and 0.50
%               If the user supplies a value of bdp smaller than 0.25, the
%               threhold found for bdp=0.25 is used.  In this case a
%               warning is produced which alerts the user that the test is
%               likely to be conservative. If on the other hand bdp is a
%               value in the interval (0.25 0.5) an average
%               between bdp=0.25 and bdp=0.5 is used (for a more refined
%               correction please see input option Tallis)
%         eff : nominal efficiency. Scalar.
%               Scalar between between 0.5 and 1-epsilon (if robest is 'MM')
%               REMARK - simulations have been done for eff = 0.85, 0.90 and
%               0.95 If the user supplies a value of eff smaller than 0.85
%               (greater than 0.95), the threshold found for eff=0.85
%               (eff=0.95) is used.  In all the other cases an average
%               is taken using the two closest values of eff
%     sizesim : simultaneous or individual size. Scalar.
%               Scalar which specifies whether simultaneous (sizesim=1) or
%               individual size is used. If sizesim is missing or equal to
%               1 a simultaneous size is used.
%     Tallis  : need to intermpolate. Scalar.
%               Scalar which has an effect just if bdp is not equal to 0.25
%               or 0.5. If Tallis=1 the program computes the ratio between
%               the asymptotic consitency factor using the breakdown point
%               supplied by the user and the closest consistency factor
%               associated to the breakdown point for which simulations
%               exist. Therefore, if for example the supplied breakdown is
%               smaller than 0.25 the program multiplies the empirical
%               threshold using bdp=0.25 by a number smaller than 1.
%               Similarly, if bdp>0.375 the program multiplies the
%               empirical threshold using bdp=0.5 by a number smaller than
%               1. If supplied bdp is very close to 0.25 or 0.5 we suggest
%               to use this option otherwise it is better to take a
%               simple average of the threholds associated to the two
%               closest breakdown points for which simulations exist. The
%               default value of Tallis is 0.
%
% Optional input arguments:
%
%
%  Output:
%
%    thresh :    Empirical threshold. Scalar.
%                Emprirical threshold which can be used in order
%                to have a test with en empirical size close to the nominal
%                size (1% individual or simultaneous)
%
%
%  More About:
%
%               We assume that the two input MAT files
%               Ind_ThreshSm.mat and Sim_ThreshSm.mat are in the same
%               folder or in the MATLAB path.
%               Ind_ThreshSm.mat contains a 3D array with the thresholds in
%               case an individual size is requested
%               Sim_ThreshSm.mat contains a
%               3D array with the thresholds in case a simultaneous size
%               is requested
%               The two 3D arrays have dimension 12-by-9-by-24
%               The 12 rows are referred to the 12 sample sizes which have
%               been considered namely n=50, 60, 70, 80, 90, 100, 150, 200,
%               250, 300, 400, 500.
%               The 9 colums are referred to the number of variables which
%               have been considered namely p=2, 3, ..., 10.
%               The third dimension is associated with the 24 estimators
%               which have been used. The order of the estimators is:
%                 ' 1'    'LTSbdp050' ;
%                 ' 2'    'LTSbdp025' ;
%                 ' 3'    'LTSrbdp050';
%                 ' 4'    'LTSrbdp025';
%                 ' 5'    'Sbdp025TB' ;
%                 ' 6'    'Sbdp050TB' ;
%                 ' 7'    'MMeff085TB';
%                 ' 8'    'MMeff090TB';
%                 ' 9'    'MMeff095TB';
%                 '10'    'Sbdp025OP' ;
%                 '11'    'Sbdp050OP' ;
%                 '12'    'MMeff085OP';
%                 '13'    'MMeff090OP';
%                 '14'    'MMeff095OP';
%                 '15'    'Sbdp025HY' ;
%                 '16'    'Sbdp050HY' ;
%                 '17'    'MMeff085HY';
%                 '18'    'MMeff090HY';
%                 '19'    'MMeff095HY';
%                 '20'    'Sbdp025HA' ;
%                 '21'    'Sbdp050HA' ;
%                 '22'    'MMeff085HA';
%                 '23'    'MMeff090HA';
%                 '24'    'MMeff095HA'.
%
%
% See also: Sreg , MMreg, LXS
%
% References:
%
%   Salini S., Cerioli A., Laurini F. and Riani M. (2014), Reliable Robust
%   Regression Diagnostics, submitted.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('robregrsize')">Link to the help function</a>
% Last modified 06-Feb-2015

%{
    % Example 1
    % Find the threshold for MM estimator, Tukey biweight rho function with
    % efficiency 0.87 (simultaneous size)
    n=232;
    p=10;
    bdp='';
    robest='MM';
    eff=0.87;
    rhofunc='TB';
    sizesim=1;
    thresh=RobRegrSize(n,p,robest,rhofunc,bdp,eff,sizesim);
%}

%{
    % Example 2
    % Find the threshold for MM estimator, take an average threhold for all
    % rho functions, and use efficiency 0.85 (simultaneous size)
    n=93;
    p=5;
    bdp='';
    eff=0.85;
    robest='MM';
    rhofunc='ST';
    sizesim=1;
    thresh=RobRegrSize(n,p,robest,rhofunc,bdp,eff,sizesim);
%}

%{
    % Example 3
    % Find the threshold for LTS estimator, use Tallis correction to infer
    % a threshold for bdp equal to 0.27 (simultaneous size)
    n=72;
    p=10;
    bdp=0.27;
    robest='LTS';
    eff='';
    rhofunc='';
    sizesim=1;
    Tallis=1;
    thresh=RobRegrSize(n,p,robest,rhofunc,bdp,eff,sizesim,Tallis);
%}

%{
    % Example 4
    % Find the threshold for S estimator and hyperbolic rho function, 
    % use Tallis correction to infer
    % a threshold for bdp equal to 0.3 (simultaneous size)
    n=100;
    p=5;
    bdp=0.3;
    robest='S';
    eff='';
    rhofunc='HY';
    sizesim=1;
    Tallis=1;
    thresh=RobRegrSize(n,p,robest,rhofunc,bdp,eff,sizesim,Tallis);
%}

%% Beginning of code

if nargin <2
        error('FSDA:RobRegrSize:missingInputs','n and/or p are missing');

else
    if p>20
        warning('FSDA:RobRegrSize:Wrongp','The number of variables is too large to produce reliable estimates')
    end
    
    if n<50 || n>500
        if n<50
            error('FSDA:RobRegrSize:Wrongn','Correction factors have been estimated for n in the range [50 500]')
        end
        
        if n>500 && p<=10
            disp('Correction factors are not needed for n greater than 500')
            error('FSDA:RobRegrSize:Wrongn','Please use the nominal chi2 threshold')
        end
        
        if n>500 && p>10
            warninging('FSDA:RobRegrSize:Wrongn','Simulations have been done for n in the range [50 500] and p<=10')
            warning('FSDA:RobRegrSize:Wrongn','The resulting interpolated factors may produce a conservative test')
            n=500;
        end
        
    end
    % Vector which contains the labels for the estimators whose size has to
    % be rescaled
    labest={'LTSbdp050','LTSbdp025','LTSrbdp050','LTSrbdp025',...
        'Sbdp025TB','Sbdp050TB','MMeff085TB','MMeff090TB','MMeff095TB',...
        'Sbdp025OP','Sbdp050OP',...
        'MMeff085OP','MMeff090OP','MMeff095OP','Sbdp025HY','Sbdp050HY','MMeff085HY',...
        'MMeff090HY','MMeff095HY','Sbdp025HA','Sbdp050HA','MMeff085HA','MMeff090HA','MMeff095HA'};
    
end


if nargin<8
    Tallis=0;
end


if nargin <3 || isempty(robest)
    robest='MM';
else
    % check if the estimator which has been provided in the list
    robestchk={'LTS','LTSr','S','MM'};
    if sum(strcmp(robest,robestchk)) ~= 1
        disp('Supplied estimator is not in the list of the estimators which have been considered')
        error('FSDA:RobRegrSize:WrongEstimator','Estimators considered are LTS, LTSr (LTS reweighted) S and MM')
    end
end

% Check
   if (strcmp(robest,'LTS') || strcmp(robest,'LTSr')) && ~isempty(rhofunc) 
       warning('FSDA:RobRegrSize:WrongEstimator','Robust estimator which has been chosen is of the LTS form')
       warning('FSDA:RobRegrSize:WrongEstimator','Therefore rho function is set to empty')
       rhofunc='';
   end     

if nargin <4 || (isempty(rhofunc) && (strcmp(robest,'S')  || strcmp(robest,'MM')))
    rhofunc='ST';
else
    if strcmp(rhofunc,'S') || strcmp(rhofunc,'MM')
        % check if rho function which has been provided is in the list
        rhofuncchk={'TB','OP','HA','HY','ST'};
        if sum(strcmp(rhofunc,rhofuncchk)) ~= 1
            disp('rho function which has been supplied is not in the list of rho functions which have been considered')
            disp('rho functions considered are')
            error('FSDA:RobRegrSize:WrongRho','TB=tukey biweight, OP=optimal, HA=Hampel, HY=hyperbolic')
        end
    end
end

if nargin<5 || isempty(bdp)
    bdp=0.5;
    corr=1;
else
    
    if bdp>0.5
        error('FSDA:RobRegrSize:WrongBdp','bdp must be a number in the interval [0+epsilon-0.5]')
    end
    
    if Tallis==1
        if bdp<0.375
            mmn=1-0.25;
            a=norminv(0.5*(1+mmn));
            T025=1-2*(1/mmn)*a*normpdf(a);
            
            mmn=1-bdp;
            a=norminv(0.5*(1+mmn));
            Tuser=1-2*(1/mmn)*a*normpdf(a);
            corr=sqrt(T025/Tuser);
        else
            mmn=1-0.5;
            a=norminv(0.5*(1+mmn));
            T05=1-2*(1/mmn)*a*normpdf(a);
            
            mmn=1-bdp;
            a=norminv(0.5*(1+mmn));
            Tuser=1-2*(1/mmn)*a*normpdf(a);
            corr=sqrt(T05/Tuser);
        end
    else
         if bdp<0.25
            disp('Simulations have been done for bdp=.25 and bdp=0.5')
            disp('The resulting interpolated factors may produce a conservative test')
            warning('FSDA:RobRegrSize:WrongBdp','Please use input option Tallis if you want to use a recalibrated threshold')
         end

        corr=1;
    end
    
end

if nargin<6 || isempty(eff)
    eff=0.90;
    if eff<0.6 || eff>=1
        error('FSDA:RobRegrSize:WrongEff','eff must be a number in the interval [0.6-(1-epsilon)]')
    end
end

% The default is to consider a simultaneous size
if nargin<7 || isempty(sizesim)
    sizesim=1;
end

% bdpchk = values of breakdown point for which simulations have been performed
bdpchk=[0.25 0.5];
% effchk = values of efficiency for which simulations have been performed
effchk=[0.85 0.9 0.95];

%% Extract from mat file the information in order to find the proper threshold
if strcmp(robest,'MM')==1
    % If efficiency is the one supplied by the user
    if intersect(eff,effchk)==eff
        % If the weight function is ST than it is necessary to take an
        % average of the threholds for all the weight function for a given
        % value of efficiency
        % posint is the boolean vector which contains the position
        % associated to the required efficiency
        if strcmp(rhofunc,'ST')
            if eff==0.85
                posint=strncmpi(labest,'MMeff085',8);
            elseif eff==0.90
                posint=strncmpi(labest,'MMeff090',8);
            elseif eff==0.95
                posint=intersect(labest,'MMeff095',8);
            end
        else
            % If the weight function is not equal to ST than it is
            % necessary to extract the position inside the MAT file
            % associated to the supplied weight function
            namefile1=[robest 'eff0' num2str(100*eff) rhofunc];
            namefile2='';
        end
    elseif  eff<0.85
        % If the used wants to have an effiency smaller than 0.85 then the
        % correction for eff=0.85 is used
        if strcmp(rhofunc,'ST')
            posint=strncmpi(labest,'MMeff085',8);
        else
            namefile1=[robest 'eff085'  rhofunc];
            namefile2='';
        end
    elseif  eff>0.85 && eff <0.9
        if strcmp(rhofunc,'ST')
            posint1=strncmpi(labest,'MMeff085',8);
            posint2=strncmpi(labest,'MMeff090',8);
            % postint = logical vector with true for the indexes in the
            % position assocaited with the MM estimators which have 0.85 or
            % 0.90 effiicency
            posint=posint1 | posint2;
        else
            % select the two efficiencies closer to those specified by the user
            namefile1=[robest 'eff085'  rhofunc];
            namefile2=[robest 'eff090'  rhofunc];
        end
        
    elseif  eff>0.90 && eff <0.95
        if strcmp(rhofunc,'ST')
            posint1=strncmpi(labest,'MMeff090',8);
            posint2=strncmpi(labest,'MMeff095',8);
            posint=posint1 | posint2;
        else
            % select the two efficiencies closer to those specified by the user
            namefile1=[robest 'eff090'  rhofunc];
            namefile2=[robest 'eff095'  rhofunc];
        end
    else
        % If the used wants to have an effiency greater than 0.95 then the
        % correction for eff=0.95 is used
        if strcmp(rhofunc,'ST')
            posint=strncmpi(labest,'MMeff095',8);
        else
            namefile1=[robest 'eff095'  rhofunc];
            namefile2='';
        end
    end
elseif strcmp(robest,'S')==1
    % If breakdown point is the one supplied by the user
    if intersect(bdp,bdpchk)==bdp
        % If the weight function is ST, it is necessary to take an
        % average of the threholds for S estimators associated to all the
        % weight functions for a given value of bdp
        % posint is the boolean vector which contains the position
        % associated to the required breakdown point inside the 3D MAT
        % file which contains all the threshold
        if strcmp(rhofunc,'ST')
            if bdp==0.25
                posint=strncmpi(labest,'Sbdp025',7);
            elseif bdp==0.5
                posint=strncmpi(labest,'Sbdp050',7);
            end
        else
            namefile1=[robest 'bdp0' num2str(100*bdp) rhofunc];
            namefile2='';
        end
        
        % if what follows is true the user has selected a bdp different from 0.25 and 0.5
    elseif  bdp<0.25 % in this case the user has selected a bdp smaller than 0.25
        if strcmp(rhofunc,'ST')
            posint=strncmpi(labest,'Sbdp025',7);
        else
            namefile1=[robest 'bdp025'  rhofunc];
            namefile2='';
        end
        
        
    else % in this case the user has selected a bdp in the interval 0.25-0.5
        if strcmp(rhofunc,'ST')
            
            if Tallis==1
                if bdp<0.375
                    posint1=strncmpi(labest,'Sbdp025',7);
                    posint2='';
                else
                    posint1=strncmpi(labest,'Sbdp050',7);
                    posint2='';
                end
            else
                posint1=strncmpi(labest,'Sbdp025',7);
                posint2=strncmpi(labest,'Sbdp050',7);
            end
            posint=posint1 | posint2;
        else
             
            if Tallis==1
                if bdp<0.375
                    namefile1=[robest 'bdp025'  rhofunc];
                    namefile2='';
                else
                    namefile1=[robest 'bdp050'  rhofunc];
                    namefile2='';
                end
            else
                % select the two breakdown points closer to those specified by the user
                namefile1=[robest 'bdp025'  rhofunc];
                namefile2=[robest 'bdp050'  rhofunc];
            end
        end
    end
elseif strcmp(robest,'LTSr')==1  % LTS reweighted
    % If breakwdon point is the one supplied by the user
    if intersect(bdp,bdpchk)==bdp
        namefile1=[robest 'bdp0' num2str(100*bdp) rhofunc];
        namefile2='';
        
    elseif  bdp<0.25
        namefile1=[robest 'bdp025'  rhofunc];
        namefile2='';
    else
        if Tallis==1
            if bdp<0.375
                namefile1=[robest 'bdp025'  rhofunc];
                namefile2='';
            else
                namefile1=[robest 'bdp050'  rhofunc];
                namefile2='';
            end
        else
            % select the two breakdown points closer to those specified by the user
            namefile1=[robest 'bdp025'  rhofunc];
            namefile2=[robest 'bdp050'  rhofunc];
        end
    end
else % RAW LTS
    % If breakdown point is the one supplied by the user
    if intersect(bdp,bdpchk)==bdp
        namefile1=[robest 'bdp0' num2str(100*bdp) rhofunc];
        namefile2='';
        
    elseif  bdp<0.25
        namefile1=[robest 'bdp025'  rhofunc];
        namefile2='';
    else
        if Tallis==1
            if bdp<0.375
                namefile1=[robest 'bdp025'  rhofunc];
                namefile2='';
            else
                namefile1=[robest 'bdp050'  rhofunc];
                namefile2='';
            end
        else
            % select the two breakdown points closer to those specified by the user
            namefile1=[robest 'bdp025'  rhofunc];
            namefile2=[robest 'bdp050'  rhofunc];
        end
    end
end

% Simulations have been performed for the values of n and p given in
% vectors nn and pp
nn=[50:10:100, 150, 200, 250, 300, 400, 500];
pp=2:10;
% rows are referred to sample size
% columns are referred to number of variables
% third dimension = estimator

% Select the slices of third dimenson of mat file which have to be used

if strcmp(rhofunc,'ST')
    % in this case Boolean vector posint has already been calculated
else
    [~,posint1]=intersect(labest,namefile1);
    [~,posint2]=intersect(labest,namefile2);
    posint=[posint1 posint2];
end

if sizesim==1
    MAT=load('Sim_ThreshSm');
    THsel=mean(MAT.ThreshAllSmsim(:,:,posint),3);
else
    MAT=load('Ind_ThreshSm');
    THsel=mean(MAT.ThreshAllSmind(:,:,posint),3);
end

if Tallis==1;
    THsel=THsel*corr;
end

THsel=[nn' THsel];
% Check is the supplied value of n is in our grid
seln=find(nn>n);
if isempty(seln); % in this case n is equal to 500
    THseln=THsel([size(THsel,1) size(THsel,1)],:);
else
    % Find the two values of n among those of the simulations which enclose
    % n. First row of THseln is associated with n smaller or equal than n
    % supplied by the user. Second row is greater than n supplied by the
    % user
    THseln=THsel(seln(1)-1:seln(1),:);
end

% If p >10 it is necessary to extrapolate by fitting a quadratic model
if p>10
    yhat=[THseln(:,1) zeros(2,1)];
    x=pp';
    x0=[1 p p.^2];
    XX=[ones(length(x),1) x x.^2];
    y=THseln(1,2:end)';
    b=XX\y;
    yhat(1,2)=x0*b;
    y=THseln(2,2:end)';
    b=XX\y;
    yhat(2,2)=x0*b;
else
    yhat=THseln(:,[1 p]);
end

% Finally interpolate to obtain the threhold for the value of n requested by the user
nU=yhat(2,1);
nL=yhat(1,1);
% Find the requested threshold by interpolation
thresh=yhat(1,2)+ ((yhat(2,2)-yhat(1,2))/(nU-nL))*(n-nL);

end
