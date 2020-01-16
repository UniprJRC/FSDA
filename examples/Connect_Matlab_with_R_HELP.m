%Connect_Matlab_with_R_HELP prepares connection with R in order to run R code from within MATLAB.
%
% REMARK: this is a third party function subject to the licence at this link:
% http://rcom.univie.ac.at/download/licenses/SC_LICENSE_HomeStudent
%
% As stated in the license, all publications based on analyses performed
% using statconnDCOM directly or indirectly (e.g. RExcel, FlexArray,
% SAM,...) must include the following citation:
%
%      Baier Thomas, & Neuwirth Erich (2007). Excel :: COM :: R.
%      Computational Statistics, Volume 22, Number 1/April 2007. Physica
%      Verlag.
%
% Based on the documentation provided by the authors of statconnDCOM and
% based on our experience, in order to connect MATLAB with R it is
% necessary to follow the next 5 steps:
%
% STEP 1: INSTALL R (>= 2.12.0)
%         When you install R in a 64bit CPU, make sure that you also
%         install the 32 bit version. In addition, click on the option
%         "Save version number in registry" during the installation.
%
% STEP 2: INSTALL R package (library) "rscproxy", required for connecting 
%         to R. Please, "TAKE SPECIAL CARE to install rscproxy into the
%         library directory of your R installation: DO NOT USE a "user"
%         library directory, if present. This has to be taken care of
%         especially on Windows Vista and Windows 7 or Windows 8. Run R (32
%         bit version, that is assuming that R3.1.1 has been installed, run
%         link Ri386 3.1.1 or run directly file Rgui.exe inside subfolder
%         i386) as Administrator. When installing package rscproxy, take
%         care to specify %R_HOME%/library as destination. For example,
%         assuming that R has been installed in C:\packages\R-3.1.2 library
%         rscproxy can be installed from R using the instruction
%         install.packages("rscproxy",lib="C:/packages/R-3.1.2/library")
%         Please, verify correct installation afterwards (open %R_HOME% in
%         Windows Explorer and verify that rscproxy has been installed into
%         %R_HOME%\library)!".%
%
% STEP 3: Download and install the Home&Student version of statconnDCOM3
%         The file can be downloaded from
%         http://rcom.univie.ac.at/download/current/statconnDCOM3.6-0B2_Noncommercial.exe
%         REMARK: you can use all the default installation options.
%
% STEP 4: Download file Testing.lic 
%         http://rcom.univie.ac.at/download/licenses/Testing.lic
%         and copy it to the installation directory of statconnDCOM
%  
% STEP 5: Download MATLAB_RLINK, to have all the matlab connection functions
%         The zip file containing all these routines, MATLAB_RLINK, can be
%         donwloaded from
%         http://www.mathworks.com/matlabcentral/fileexchange/5051-matlab-r-link.
%         Once the .zip file has been downloaded, extract it to a folder
%         and add such folder to the MATLAB path. For example if the
%         content of MATLAB_R_LINK has been extracted to
%         D:\matlab\MATLAB_R_LINK, to add this folder to the path it is
%         necessary to run the following line
%         addpath('D:\matlab\MATLAB_R_LINK'); To add this folder
%         permanently to the path click Home|Set Path and then click on
%         Save, or alternatively from the prompt write savepath.
%
%         IN ORDER TO TEST THAT EVERYTHING IS OK IT, IS POSSIBLE TO USE THE
%         THREE EXAMPLES BELOW
%
%         REMARK: if when you run openR from the MATLAB command prompt the
%         following line appears
%
%         Error using openR (line 68)
%         Cannot connect to R.
%         Error using COM.StatConnectorSrv_StatConnector/Init
%         Error: Object returned error code: 0x80040013
%
%         there can be one of the following explanations:
%         1) package rscproxy has not been installed under %R_HOME%\library
%         2) R has been installed without the option "Save version number
%            in registry". In this last case it is necessary to reinstall R
%
%         REMARK: if just example 1 below can be executed but not example 2
%         and the result of instruction rn1 = evalR('rnorm(1)') is an
%         empty value, or alternatively if evalR('rnorm(1)') produces the
%         following error
%
%         "Invoke Error, Dispatch Exception: Object is static; operation
%         not allowed"
%
%         this is due to the fact that the connection software cannot find
%         R dll library (called Rlapack.dll which is generally installed in
%         %R_Home%\bin\i386.
%
%         In some computers it may be due to the fact that R had been
%         previously installed in folder "C:\Program files". In this case it
%         is necessary to reinstall R in a different folder without spaces.
%
%         Alternatively, it is necessary to edit file SC_Proxy.cfg located
%         in statconn\DCOM and add the link to file Rlapack.dll
%
%         If the problem persists it is necessary to install statconnDCOM
%         using administrator privileges in a path different from
%         "C:\Program files" and not to use spaces in the path.
%
%         An additional test which can be done is from Start|All
%         Programs|statconn|DCOM|Server 01 - Basic Test 
%         If, after clicking on the R button in the Gui, the following
%         message appears:
%
%          Code: -2147221455
%          Text: licensing: license expired
%          Releasing StatConnector Server...Done
%   
%         it means that  file Testing.lic (see STEP 2 above) has not been
%         copied
%
% Copyright 2008-2019.
% FSDA toolbox
%
%
%$LastChangedDate::                      $: Date of the last commit

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

%% EXAMPLE 2: CHECK THE GENERATION OF RANDOM NUMBERS FROM R

openR

% ask R what is the type of random number generator
a = evalR('ok <- RNGkind()');
% set the seed based on the Mersenne-Twister generator
evalR('set.seed(1234, kind=''Mersenne-Twister'' , normal.kind = ''Inversion'')');
% the random streem
b = evalR('.Random.seed');
% generate a random number
rn1 = evalR('rnorm(1)');
disp(rn1)
% This is what you should get from the MATLAB prompt
% >>rn1
% 
% rn1 = -1.2071

% set  again the seed
evalR('set.seed(1234, kind=''Mersenne-Twister'' , normal.kind = ''Inversion'')');
% generate the same random number as before
rn2 = evalR('rnorm(1)');
% >>rn2
% 
% rn2 = -1.2071

% a different random number
rn3 = evalR('rnorm(1)');
% >>rn3
% 
% rn3 = 0.2774

closeR

%% EXAMPLE 3: CHECK A R DEMO

openR
evalR('demo("persp")');
closeR

% Other examples can be found in files Rdemo.m and Rdemo.html in folder MATLAB_R_LINK
% At the MATLAB prompt type
% which('Rdemo') or which('Rdemo.html')
% to see where these file are located
