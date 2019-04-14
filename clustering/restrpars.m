function [Sigma]  = restrpars(SigmaB, niini, pa)
% SigmaB = v-times-v-times-k = empirical covariance matrix
% pa = structure containing modeltype, number of iterations .....
Sigma=SigmaB;
K=length(niini);
lmd=NaN(1,K);

% OMG = initialize 3D array containing rotation matrices
OMG=zeros(size(Sigma));

% Tolerance associated to the maximum of the elements which have to be
% constrained. If the maximum  is smaller than zerotol restreigen
% procedures returns in output what has been given in input. For example,
% if the elements which have to be constrained are the eigenvalues of the
% covariance matrices and the max of the eigenvalues is smaller than
% zerotol it means that all n points are concentrated in k points and there
% is a perfect fit therefore no further changes on the eigenvalues is
% required.
zero_tol=pa.zerotol;

% pa.pars = character vector with three letters specifying the type of the
% 14 contraints (i.e. EEE, CVVV, EVE, ...)
pars=pa.pars;


% VEV VEE EVE VVE require iterations.
% All the other specification do not
if strcmp(pars,'EEE') || strcmp(pars,'VVV') || strcmp(pars,'EVV') ...
        || strcmp(pars,'EEV') || strcmp(pars(3),'I')
    pa.iterDSR = 1;
end

% If equal determinants are imposed pa.cdet=1
if strcmp(pars(1),'E')
    pa.cdet = 1;
end

% If OMG is identity shape restriction parameter within gouprs is set to 1
if strcmp(pars(2),'I')
    pa.shw = 1;
end

%% Initialization part
if strcmp(pars(3),'E')
    % In the common principal components case it is necessary to find
    % initial values for OMG, GAM and lmd
    
    % Find initial values of lmd (unconstrained determinants) and OMG (rotation)
    [lmd, OMG]  = initR(SigmaB, niini, pa);
    % Find initial constrained shape matrix GAM 
    % pa.shw and pa.shb constraining parameters are used
    [GAM] =restrshapepars(lmd, OMG, SigmaB, niini, pa);
    % Find initial constrained determinants (lmd vector)
    lmd =restrdeterpars(GAM, OMG, SigmaB, niini, pa, zero_tol);
    
elseif  strcmp(pars(3),'V')
    % Find initial (and final value for OMG)
    for k=1:K
        [V,~]= eig(SigmaB(:,:,k));
        V=fliplr(V);
        % D=fliplr(D);
        OMG(:,:,k)=V;
    end
    
else % The remaning case is when **I
    % Find initial (and final value for OMG). 
    % in this case OMG is a 3D arry contaning identity matrices
    eyep=eye(pa.p);
    for k=1:K
        OMG(:,:,k)=eyep;
    end
end

% End of initialization

%% Beginning of iterative process
for i=1:pa.maxiterDSR
    
    % In the **E case it is necessary to update in each step of the
    % iterative procedure OMG
    % Variable shape: update OMG (rotation)
    if strcmp(pars,'VVE') || strcmp(pars,'EVE')
        Omega_=OMG(:,:,1);
        % parameter iterR is used here
        [OMG]  = cpcV(lmd, GAM, Omega_, SigmaB, niini, pa);
    elseif strcmp(pars,'VEE') || strcmp(pars,'EEE')
        % Equal shape: update OMG (rotation)
        % old eiggiroe
        [OMG]  = cpcE(lmd, SigmaB, niini, pa);
    else
        % In all the other cases OMG is not updated
    end
    
    % Find new value of GAM (shape matrix)
    [GAM] =restrshapepars(lmd, OMG, SigmaB, niini, pa);
    % Find new value of lmd (determinants)
    [lmd]=restrdeterpars(GAM, OMG, SigmaB, niini, pa, zero_tol);
end

% Check if all is well
codeZero = max(max(GAM)) > zero_tol;

if  ~codeZero
    return
end

% Reconstruct the cov matrices using final values of lmd, OMG and GAM
for k=1:K
    Sigma(:,:,k) = lmd(k) * OMG(:,:,k)* diag(GAM(:,k))* (OMG(:,:,k)');
end

end