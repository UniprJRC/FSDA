%% PREPARE CONNECTION WITH R in order to execute R commands from within MATLAB.
%
%   In order to achieve the goal above it is necessary to follow the 3
%   following steps.
%
% STEP 1: INSTALL R-(D)COM Interface
%         The installation file is called statconnDCOM.latest.exe and can
%         be downloaded from
%         http://rcom.univie.ac.at/download/current/statconnDCOM.latest.exe
%         Remark: statconnDCOM.latest.exe is the updated version of the
%         installation program which can be downloaded
%         from http://lib.stat.cmu.edu/R/CRAN/other-software.html (which
%         dows not seem to work)
%         REMARK: when you install statconnDCOM.latest.exe setup checks
%         whether there are previous installations. If this is the case it is
%         necessary to remove them.
%
% STEP 2: INSTALL R package (library) "rscproxy", required for connecting to R
%         As explained in the setup of the installation statconnDCOM.latest.exe
%         "TAKE SPECIAL CARE to install rscproxy into the library directory
%         in your R installation and not in some "user" library directory.
%         This has to be taken care of especially on Windows Vista and
%         Windows 7 or Windows 8. Run R as Administrator when installing
%         these packages and take care to specify %R_HOME%/library as
%         destination. Please verify correct installation afterwards (open
%         %R_HOME% in Windows Explorer and verify that rscproxy has been
%         installed into %R_HOME%\library)!".
%
% STEP 3: Download MATLAB_RLINK, to have all the matlab connection functions
%         The zip file contanining all these routines is called
%         MATLAB_RLINK and can be donwloaded from the link
%         http://www.mathworks.com/matlabcentral/fileexchange/5051-matlab-r-link.
%         Once the .zip file has been downloaded extract it to a folder and
%         add to the MATLAB path the folder where you extracted the .zip
%         file.
%         For example if the content of MATLAB_R_LINK has been extracted to
%         C:\matlab\MATLAB_R_LINK, to add this folder to the path it is
%         necessary to run the following line
%         addpath('C:\matlab\MATLAB_R_LINK');
%         To add this folder permanently to the path from Home|Set PAth and
%         clicj on Save
%
%         IN ORDER TO TEST THAT EVERYTHING IS OK from
%               Start | All programs | statconn | DCom
%         it is possible to run
%               Server 01 - Basic Test
%
%         Alternatively run openR from the MATLAB command prompt
%
%         REMARK: if the following line appears
%         Error using openR (line 68)
%         Cannot connect to R.
%         Error using COM.StatConnectorSrv_StatConnector/Init
%         Error: Object returned error code: 0x80040013
%
%         there can be one of the following explanations
%         1) R_HOME is not the windows path. If this is the case it is
%         necessary to add it to the path.
%         2) package rscproxy has not been installed under the %R_HOME%\library
%            folder
%         3) R has been installed without the option
%          Save version number in registry". In this last case it is necessary
%          to reinstall R
%

%% EXAMPLE 1: create from within R an identity matrix

openR
%
% run a R command, e.g. the generation of a random number
rn = evalR('diag(3)');
%
% This is what you should get
%
%         rn =
%
%              1     0     0
%              0     1     0
%              0     0     1
%
%
% close the R connection
closeR %close the connection

%% EZAMPLE 2: CHECK THE GENERATION OF RANDOM NUMBERS FROM R

openR

% ask R what is the type of random number generator
a = evalR('ok <- RNGkind()');
% set the seed based on the Mersenne-Twister generator
evalR('set.seed(1234, kind=''Mersenne-Twister'' , normal.kind = ''Inversion'')');
% the random streem
b = evalR('.Random.seed');
% generate a random number
rn1 = evalR('rnorm(1)');
% This is what you should get from the MATLAB prompt
% >>rn1
% 
% rn1 =
% 
%    -1.2071
% set  again the seed
evalR('set.seed(1234, kind=''Mersenne-Twister'' , normal.kind = ''Inversion'')');%Super
% generate the same random number as before
rn2 = evalR('rnorm(1)');
% >>rn1
% 
% rn2 =
% 
%    -1.2071
% a different random number
rn3 = evalR('rnorm(1)');
% >>rn3
% 
% rn3 =
% 
%     0.2774
closeR


