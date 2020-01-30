function [ypred, yci] = predci(X,beta,Sigma,mse,dfe,alpha,sim,pred,hasintercept)
%preci computes prediction intervals in the linear regression model
%
%  Required input arguments:
%
%    X:         Points at which routine predci predicts responses. Matrix. 
%               The size of X is k x p where k is the number of points for
%               which confidence interval is requested and p is the number
%               of parameters in the regression model
%    beta:      Vector of regression coefficients. Column vector. Column
%               vector containing the estimated regression coefficients
%    Sigma:     Covariance matrix of coefficient estimates. Matrix.
%               Sigma is equal to mse*inv(Xori'*Xori) where Xori is the
%               matrix of explanatory variables which was used to obtain
%               beta
%     mse :     Mean squared error (residuals). Scalar. mse  = SSE/dfe,
%               where SSE is the sum of squares of residuals
%     dfe :     Degrees of freedom for error (residuals). Scalar.
%               dfe is equal to the number of observations (number of rows
%               of matrix Xori) minus the number of estimated coefficients.
%    alpha:     Positive scalar from 0 to 1. Confidence level of yci is 100(1 – alpha)%.
%               For example 0.05 implies a 95% confidence interval.
%    sim  :     Reference distribution. Logical value.
%               Logical value specifying whether the confidence bounds are
%               for all predictor values simultaneously (true), or hold for
%               each individual predictor value (false). Simultaneous
%               bounds are wider than separate bounds, because it is more
%               stringent to require that the entire curve be within the
%               bounds than to require that the curve at a single predictor
%               value be within the bounds. if sim is true confidence bands
%               are based on the F distribution else they are based on the
%               t distribution
%    pred :     Boolean specifying the type of prediction. Boolean. If pred
%               is true it computes prediction interval for new
%               observations. This results in wider bounds because the
%               error in a new observation is equal to the error in the
%               estimated mean value, plus the variability in the
%               observation from the true mean else if pred is false it
%               produces  confidence bounds for the fitted mean values.
% hasintercept: Boolean associated to intercept. Boolean. 
%              Option hasintercept is just used only if sim and pred are
%              both true
%
%
%
% Output
%
% ypred = Compute the predicted values at the new X.
% yci = Confidence intervals. A two-column matrix with each row providing
% one interval. The meaning of the confidence interval depends on the
% settings of input parameters sim and pred
%

% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

ypred = X * beta;

if nargout > 1 % Calculate confidence interval
    
    if (pred) % prediction interval for new observations
        varpred = sum((X*Sigma) .* X,2) + mse;
    else % confi interval for fitted curve
        varpred = sum((X*Sigma) .* X,2);
    end
    
    if (sim) % simultaneous
        if (pred)
            % For new observations.
            if (hasintercept)
                % Jacobian has constant column.
                sch = length(beta);
            else
                % Need to use a conservative setting.
                sch = length(beta) + 1;
            end
        else
            % For fitted curve.
            sch = length(beta);
        end
        crit = sqrt(sch * finv(1-alpha, sch, dfe));
    else % pointwise
        crit = tinv(1-alpha/2,dfe);
    end
    delta = sqrt(varpred) * crit;
    yci = [ypred-delta ypred+delta];
end
end

