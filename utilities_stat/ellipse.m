function ellipse(mu,Sigma,conflev)
%ellipse generates an ellipse given mu (location vector) and Sigma (scatter matrix)
%
%   The ellipse is generated using the equation:
%
%    (x-mu)' \Sigma^{-1} (x-mu) = c
%
%<a href="matlab: docsearchFS('ellipse')">Link to the help function</a>
%
% Required input arguments:
%
% mu    : vector of two elements associated with the center of the ellipse
% Sigma : 2 x 2 symmetric matrix 
%
% Optional input arguments:
%
%
%               conflev : scaalar, confidence level (if c is not specified the value 
%                   chi2inv(0.95,2) is used
%               h : the axis handle of a figure where to send the ellipse.
%                   This can be used to host the ellipse in a subplot of
%
% See also ellipsoid
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ellipse')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    rho=-2;
    A=[4 rho; rho 3 ];
    mu=[1.5 1];
    ellipse(mu,A);
%}

%% Beginning of code
% If the user has provided has input a column vector take the transpose
if ~isrow(mu)
    mu=mu';
end

if nargin<3;
    c = chi2inv(0.95,2);
else
    c = chi2inv(conflev,2);
end

% Compute eigenvalues and eigenvectors of matrix Sigma
% Set to 0 elments smaller than 1e-14 to avoid numerical problems with
% computation of eigenvalues
Sigma(abs(Sigma(:))<1e-14)=0;

[Gam,Lam] = eig(Sigma*c);

% Make sure that Lam(1,1) is smaller than Lam(2,2);
if Lam(1,1)>Lam(2,2);
    Gam=Gam(:,[2 1]);
    Lam1=zeros(2,2);
    Lam1(1,1)=Lam(2,2);
    Lam1(2,2)=Lam(1,1);
    Lam=Lam1;
end

%disp(Gam)
 for j=1:2;
     if Gam(1,j)<0 && Gam(2,j)<0;
        Gam(:,j)=-Gam(:,j);
     else   
        if Gam(1,j)*Gam(2,j)<0;
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
Xori=bsxfun(@plus,X, mu);

% hold('on')
plot(Xori(:,1),Xori(:,2),'r');

% Add line associated with major axis
ax1=[-lenax1 0; lenax1 0];
ax1ori=ax1*Gam;
ax1ori=bsxfun(@plus,ax1ori, mu);
line(ax1ori(:,1),ax1ori(:,2),'Color','r');

% Add line associated with minor axis
ax2=[0 -lenax2;0  lenax2];
ax2ori=ax2*Gam;
ax2ori=bsxfun(@plus,ax2ori, mu);
line(ax2ori(:,1),ax2ori(:,2),'Color','r');

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
