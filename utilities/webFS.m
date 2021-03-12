function webFS(varargin)
% In case of MATLAB 2021A if the user has chosen as Documentation Location
% "web on Mathworks" and the user is searching for a FSDA function temporarily
% change the Domain in order to prevent to search in Mathworks web site.
com.mathworks.mlservices.MLHelpServices.setDocCenterDomain('http://rosa.unipr.it');
web(varargin{:});
com.mathworks.mlservices.MLHelpServices.setDocCenterDomain('https://www.mathworks.com/');
end
