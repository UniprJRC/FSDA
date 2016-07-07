function sigmaest = eigs_sigma(y,varargin)
%eigs_sigma sets a non zero value for the optional parameter sigma of function eigs
% Function eigs, when sigma is omitted, returns the eigenvalues largest in
% magnitude; when sigma is zero, the eigenvalues smallest in magnitude are
% found. However, when there are duplicated rows in the input matrix,
% option sigma=0 produces an error which can be circumvented by changing
% zero into an appropriately small positive value. Given input data y,
% function eigs_sigma choses such a positive sigma value.
%
%
% See also: eigs
%
% Copyright 2008-2016.
% Written by FSDA team
%
% Last modified 14-06-2016
%
% Examples
%
%{

    clear all;
    % a symmetric matrix
    A = gallery('moler',5,0.3);
    disp('larger eigenvalues');
    lambda = eigs(A)
    disp('smaller eigenvalues');
    sigma = 0;
    lambda = eigs(A,size(A,1),sigma)

%}

%{ 
    % example_producing_error.
    % the same matrix, but with two identical raws
    clear all;
    A = gallery('moler',5,0.3);
    A(5,:)=A(4,:);
    disp('larger eigenvalues are still fine');
    lambda = eigs(A)
    disp('but smaller eigenvalues produce an error');
    sigma = 0;
    lambda = eigs(A,size(A,1),sigma)
%}

%{ 

% the same matrix, but with two identical raws
clear all;
A = gallery('moler',5,0.3);
A(5,:)=A(4,:);
disp('larger eigenvalues are still fine');
lambda = eigs(A)
disp('using eigs_sigma smaller eigenvalues do not produce an error');
sigma = eigs_sigma(A);
lambda = eigs(A,size(A,1),sigma)

%}


%% Beginning of code
if nargin<2
    msg = 0;
else
    msg =1;
end

% define different levels of tolerance for convergence
smalltol = eps('double');
bigtol = eps('single');

% lu factorization: [L,U,pp] = lu(y,'vector') returns an upper triangular
% matrix in U, a lower triangular matrix L with a unit diagonal, and a
% permutation vector pp, such that y(pp,:)=L*U.
[L,U,~] = lu(y,'vector');
dU = diag(U);

sigmaest = 0;
if any(dU == 0) || any(diag(L) == 0)
    sigmaest = smalltol;
    if msg
        shiftSingular = 'Sigma=0 is an exact eigenvalue: \nI try a bigger value for sigma.\n\n';
        fprintf(shiftSingular);
    end
end
if any(dU ~= 0) && (min(abs(dU)) / max(abs(dU)) < smalltol)
    
    sigmaest = bigtol;
    
    if msg
        shiftNearSingular = 'Sigma is near an exact eigenvalue: \nThe algorithm may not converge, thus I try a bigger value for sigma.\n\n';
        fprintf(shiftNearSingular);
    end
end

end
