function [Ell , he] = ellipse(mu, Sigma, conflev, Color, axesellipse)
%ellipse generates an ellipse given mu (location vector) and Sigma (scatter matrix)
%
%
%<a href="matlab: docsearchFS('ellipse')">Link to the help function</a>
%
%   The ellipse is generated using the equation:
%   \[
%    (x-\mu)' \Sigma^{-1} (x-\mu) = c_{conflev}^2.
%   \]
%    The length of the i-th principal semiaxis $(i=1, 2)$ is $c \lambda_i$ where
%    $\lambda_i$ is an eigenvalue of $\Sigma$.
%
% Required input arguments:
%
% mu    : Center of the ellipse. Vector.
%         Vector with two elements associated with the center of the
%         ellipse
% Sigma : 2 x 2 symmetric positive definite matrix. Matrix. Inverse of the
%         matrix of the quadratic form which defines the equation of the
%         ellipse. Sigma is interpretable as the covariance matrix of the
%         original data points.
%
% Optional input arguments:
%
%       conflev : Confidence level. Scalar.
%                 Confidence level which controls the size of the ellipse.
%                 If conflev is not specified the value
%                 chi2inv(0.95,2) is used.
%                    Example - 'conflev', 0.99
%                    Data Types - single | double
%
%       Color   : LineColor of the ellipse. String or 3 elements numeric vector.
%                 Line color, specified as an RGB triplet, a color
%                 string, or 'none'. If you specify the Color as
%                 'none', then the line is invisible.
%                 An RGB triplet is a three-element row vector whose
%                 elements specify the intensities of the red,
%                 green, and blue components of the color. The
%                 intensities must be in the range [0,1], for
%                 example, [0.4 0.6 0.7].
%                    Example - 'Color', 'r'
%                    Data Types - [0 0 1] (default) | RGB triplet | color string | 'none'
% axesellipse   : axes of the ellipse. Boolean. If axes is true (default)
%                 dottted lines along the major axes of the ellipse are
%                 drawn else just the ellipse contour appears.
%                    Example - 'axesellipse', false
%                    Data Types - Boolean
%
%
%  Output:
%
%       Ell   :   x and y coordinates of the ellipse. Matrix.
%                 630-by-2 matrix containing the x and y coordinate of the
%                 ellipse.
%                 1st column = x coordinates;
%                 2nd column = y coordinates.
%       he   :    Vector of chart line objects. matlab.graphics.chart.primitive.Line.
%                 A column vector of chart line objects. It can be used to
%                 modify properties of a specific chart line of the plot
%                 containing the ellipse after it is created. For a list of
%                 properties, see Chart Line Properties.
%
%
% See also: ellipsoid
%
% References:
%
%   Mardia, K.V., J.T. Kent, J.M. Bibby (1979). Multivariate Analysis. Academic
%   Press, London-New York-Toronto-Sydney-San Francisco.
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ellipse')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Example using all default options.
    rho=-2;
    A=[4 rho; rho 3 ];
    mu=[1.5 1];
    ellipse(mu,A);
%}

%{
    %% Draw the ellipse using a blue color line.
    close all
    rho=-2;
    A=[4 rho; rho 3 ];
    mu=[1.5 1];
    Color=[0 0 1];
    ellipse(mu,A,[],Color);
%}

%{
    %% Draw an ellipse and fill it with yellow color.
    close all
    rho=-2;
    A=[4 rho; rho 3 ];
    mu=[1.5 1];
    Ell=ellipse(mu,A);
    patch(Ell(:,1),Ell(:,2),'y');
%}

%{
    %% 99 per cent confidence ellipse.
    % Generate 1000 bivariate normal data and add the ellipse which
    % contains approximately 990 of them.
    rng('default')
    rng(20)     % For reproducibility
    % Define mu and Sigma
    mu = [2,3];
    Sigma = [1,1.5;1.5,3];
    Y = mvnrnd(mu,Sigma,1000);
    figure
    hold on;
    plot(Y(:,1),Y(:,2),'o');
    % add an ellipse to these points
    Ell=ellipse(mu,Sigma,0.99);
    axis equal
    % Count number of points inside the ellipse
    disp('Number of points inside the ellipse')
    disp(sum(inpolygon(Y(:,1),Y(:,2),Ell(:,1),Ell(:,2))))
%}

%{
    % 99 per cent confidence ellipse without showing the major axes.
    % First close all previous plots
    close all
    % Generate 1000 bivariate normal data and add the ellipse which
    % contains approximately 990 of them.
    rng('default')
    rng(20)     % For reproducibility
    % Define mu and Sigma
    mu = [2,3];
    Sigma = [1,1.5;1.5,3];
    Y = mvnrnd(mu,Sigma,1000);
    figure
    hold on;
    plot(Y(:,1),Y(:,2),'o');
    % add an ellipse to these points and do not show the major axes
    AxesEllipse = false;
    Ell=ellipse(mu,Sigma,0.99,[],AxesEllipse);
    axis equal
    % Count number of points inside the ellipse
    disp('Number of points inside the ellipse')
    disp(sum(inpolygon(Y(:,1),Y(:,2),Ell(:,1),Ell(:,2))))
%}



%% Beginning of code
% Specify line width of the ellipse and of its axes
LineWidth=2;

% If the user has provided has input a column vector take the transpose
if ~isrow(mu)
    mu=mu';
end

if nargin<3 || isempty(conflev)
    c = chi2inv(0.95,2);
else
    c = chi2inv(conflev,2);
end

% Use default black color
if nargin<4 || isempty(Color)
    Color = [ 0 0 0];
end


% Compute eigenvalues and eigenvectors of matrix Sigma
% Set to 0 elments smaller than 1e-14 to avoid numerical problems with
% computation of eigenvalues
Sigma(abs(Sigma(:))<1e-14)=0;

[Gam,Lam] = eig(Sigma*c);

% Make sure that Lam(1,1) is smaller than Lam(2,2);
if Lam(1,1)>Lam(2,2)
    Gam=Gam(:,[2 1]);
    Lam1=zeros(2,2);
    Lam1(1,1)=Lam(2,2);
    Lam1(2,2)=Lam(1,1);
    Lam=Lam1;
end

%disp(Gam)
for j=1:2
    if Gam(1,j)<0 && Gam(2,j)<0
        Gam(:,j)=-Gam(:,j);
    else
        if Gam(1,j)*Gam(2,j)<0
            if Gam(1,j)*Sigma(1,2)>0
                Gam(:,j)=-Gam(:,j);
            end
        end
    end
end
% Note that the eigenvalues of matrix Sigma^-1
% are simply 1/Lam(1,1) and 1/Lam(2,2)
% The length of the semiaxis of the ellipse are simply
%  sqrt(Lam(1,1)) and sqrt(Lam(2,2));

th=(0:0.01:(2*pi+0.01))';
%plot(3*cos(th),sin(th));
% lenax1=(1/sqrt(Lam(1,1)));
% lenax2=(1/sqrt(Lam(2,2)));
lenax1=sqrt(Lam(1,1));
lenax2=sqrt(Lam(2,2));
%lenax1=sqrt(Lam(2,2));
%lenax2=sqrt(Lam(1,1));


xx=lenax1*sin(th);
yy=lenax2*cos(th);
X=[xx yy]*Gam;

% Add the means
Ell=bsxfun(@plus,X, mu);

% hold('on')
he = plot(Ell(:,1),Ell(:,2),'Color',Color,'LineWidth',LineWidth);

if nargin<5 || (nargin ==5 &&   axesellipse == true)
    
    % Add line associated with major axis
    ax1=[-lenax1 0; lenax1 0];
    ax1ori=ax1*Gam;
    ax1ori=bsxfun(@plus,ax1ori, mu);
    line(ax1ori(:,1),ax1ori(:,2),'Color',Color,'LineWidth',LineWidth-1,'LineStyle','--');
    
    % Add line associated with minor axis
    ax2=[0 -lenax2;0  lenax2];
    ax2ori=ax2*Gam;
    ax2ori=bsxfun(@plus,ax2ori, mu);
    line(ax2ori(:,1),ax2ori(:,2),'Color',Color,'LineWidth',LineWidth-1,'LineStyle','--');
end

% axis equal


% disp(Gam);
% disp(Lam);
% disp(Sigma)
% disp('%%%%%%%%%%%%%%%%%%%')
%disp(Gam)

%% Alternative code

% deter=Sigma(1,1)*Sigma(2,2)-Sigma(1,2)^2;
%
% ylimit=sqrt(chi2inv(conflev,2)*Sigma(2,2));
%
% y=-ylimit:0.005*ylimit:ylimit;
% sqtdi=sqrt(deter*(ylimit^2-y.^2))/Sigma(2,2);
% sqtdi([1,end])=0;
% b=mu(1)+Sigma(1,2)/Sigma(2,2)*y;
% x1=b-sqtdi;
% x2=b+sqtdi;
% y=mu(2)+y;
% coord=[x1,x2([end:-1:1]);y,y([end:-1:1])]';
% plot(coord(:,1),coord(:,2))

end
%FScategory:UTISTAT
