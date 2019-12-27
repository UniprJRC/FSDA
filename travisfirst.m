% add appropriate path
% run addFSDA2path.m  
% LXS with default input and output.
    % Compute LMS estimator without reweighting, add intercept to matrix X
    % and do not produce plots.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)+6;
  %   [out]=LXS(y,X);
