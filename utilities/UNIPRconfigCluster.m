function UNIPRconfigCluster(user, varargin)
%UNIPRconfigCluster creates a cluster profile for MATLAB Parallel Server
%
%
%<a href="matlab: docsearchFS('UNIPRconfigCluster')">Link to the help function</a>
%
%
% The routine is designed to use the High Performance Computing Facilities
% of the University of Parma (HPCUNIPR),
% https://www.hpc.unipr.it/dokuwiki/doku.php?id=calcoloscientifico:userguide
% however, it is included in the general utilities of FSDA, because it can
% be easily extended to exploit other cloud services. UNIPRconfigCluster
% creates a cluster profile for MATLAB Parallel Server on HPCUNIPR and
% enables the current client (Windows) PC to become a member of the remote
% cluster.
%
% PREREQUISITES: this routine assumes that: 
% 1) The user has an account on the
%   remote server (in this case the HPCUNIPR), has generated a SSH key pairs
%   and has appended his/her public key on the remote server inside the file
%   named "authorized_keys" which is located in the .ssh
%   folder of the user. For more information to have/request an account on the
%   HPCUNIPR cluster, please see the userguide at the web address
%   https://www.hpc.unipr.it/dokuwiki/doku.php?id=calcoloscientifico:userguide
% 2) MATLAB addon "Parallel Computing
%   Toolbox plugin for MATLAB Parallel Server with Slurm" has been
%   installed on the local computer and that the parallel computing toolbox
%   is also installed.
% 3) In the local computer you have mapped a network drive to
%   your directory in the remote system. For example in the case of the
%   University of Parma the remote directory (assuming the user is
%   paolo.andrei) is \\sshfs\paolo.andrei@gui.hpc.unipr.it. On windows
%   systems the software which enable to do this mapping are called
%   SSHFS-Win · SSHFS for Windows and can be downloaded from the github
%   address https://github.com/billziss-gh/sshfs-win. The two .msi files to
%   install are called sshfs-win-3.5.20357-x64.msi and
%   sshfs-win-3.5.20357-x86.msi.
%   Note that the default letter of the network drive is H: but it can be
%   changed using option MapNetworkDrive
% 4) In the remote system you have created the path ~/parallel/matlab
% 
% Using the varargin it is possible to specify the typical parameters
% of the SLURM workload manager (https://slurm.schedmd.com/), that is the
% requested number of workers, number of threads for CPU, number of tasks
% per node, the amount of RAM ....
%
% Required input arguments:
%
%    user:      user SSH credentials on the remote cluster. Character or string.
%               User must be specified as a char or as a string. In the
%               case of UNIPR user must be in the format of name.surname
%               (same format of the e-mail account without the domain
%               part). 
%               Example - 'paolo.andrei' or "paolo.andrei"
%               Data Types - character or string
%
% Optional input arguments:
%
%
%  NumWorkers :  number of workers per CPU. Positive integer.
%               Number of cores to use to perform the calculations.
%               The default value is 2.
%               Example - 'NumWorkers', 2
%               Data Types - double
%
% NumThreads:   number of threads per CPU. Positive integer. 
%               The number of threads that each core is going to use. The
%               default value is 2. 
%               Example - 'NumThreads', 2
%               Data Types - double
%
%   tasksxnode: number of tasks per node. Positive integer. 
%               The default value is 1. Note that MATLAB is able to use
%               multiple CPU-cores via libraries that have been written
%               using shared-memory parallel programming models like
%               OpenMP, Intel Threading Building Blocks (TBB) or pthreads.
%               For pure multithreaded codes, only a single node and single
%               task can be used and the optimal value of cpus-per-task is
%               sought. For more information on tasks per node please see
%               https://researchcomputing.princeton.edu/support/knowledge-base/scaling-analysis
%               Example - 'tasksxnode', 2
%               Data Types - double
%
%   memxcpu:    amount of RAM per CPU expressed in GigaBytes. 
%               Character or string containing a positive integer followed by G.
%               The default value is '2G'.
%               Note that once you are logged into the cluster it is
%               possible to use the command sinfo in order to find maximum
%               CPU/memory per node. 
%               https://support.ceci-hpc.be/doc/_contents/SubmittingJobs/SlurmFAQ.html
%               Example - 'memxcpu','2G' or 'memxcpu',"8G"
%               Data Types - char or string
%
% IPaddress: initialize and connect the Windows client to the HPC cluster. 
%               Logical, or char/string containing the IPaddress which must
%               be added to the remote cluster. 
%               The default value is true, that is the routine uses the
%               IPaddress of the local windows client. If IPaddress
%               is false only the SLURM parameters are set and in order to
%               to connect the Windows client to the cluster it is necessary to run
%               routine UNIPRruncluster. Note that routine
%               UNIPRConfigCluster can be launched every time it is
%               necessary to change the slurm parameters. If the SLURM
%               parameters do not change, than it is enough to call routine
%               UNIPRruncluster to add the current local machine to the
%               cluster.
%               Example - 'IPaddress', false
%               Data Types - boolean or character or string containing a valid IP address
%
%   MapNetworkDrive : letter which identifies network drive. Character or string.
%               Letter followed by : which identifies the mapped network
%               drive to your directory in the remote system. The default
%               value of MapNetworkDrive is 'H:'.
%               Example - 'MapNetworkDrive', 'Z:'
%               Data Types - character or string
%
%  Output:
%
% See also:     UNIPRruncluster
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('UNIPRconfigCluster')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%

%{
    % Use all default options.
    % We assume that john.doe is a valid user inside UNIPR HPC.
    UNIPRconfigCluster('john.doe')
%}

%{
    % Use varargin to specify SLURM paramters
    UNIPRconfigCluster('john.doe','tasksxnode', 1,...
    'memxcpu','2G', 'NumWorkers', 2, 'NumThreads', 2)
%}

%{
    % Create or update the MATLAB parallel profile only.
    % In this example just the MATLAB cluster profile is created and
    % routine UNIPRrunCluster is not called.
    UNIPRconfigCluster('john.doe', 'IPaddress', false)
%}


%{
    %  Create or update the MATLAB parallel profile and add local client to the remote cluster.
    UNIPRconfigCluster('john.doe', 'IPaddress', true)
    % Remark: given that the default option IPaddess is set to true
    % the previous instruction was equivalent to the instruction
    % UNIPRconfigCluster('john.doe')
%}

%{
    % run also the function runclusterUNIPR() with specific IP address.
    % With option IPaddress it is possible to control the IP address with
    % which the user wants to join the remote cluster.
    UNIPRconfigCluster('john.doe', 'IPaddress', '192.168.1.10')
%}

%% Beginning of code

if verLessThanFS([9 11])
    error('FSDA:ConfigClusterUNIPR:WrongMATLABVersion','At least MATLAB R2021b is needed.');
end

if nargin < 1
    error('FSDA:ConfigClusterUNIPR:WrongInputOpt','user name in the format name.surname is missing.');
end

% make sure that user is in char format
user=char(user);

if ~contains(user, '.')
    error('FSDA:ConfigClusterUNIPR:WrongInput','user name must be in the format name.surname.');
end

if contains(user, '@')
    error('FSDA:ConfigClusterUNIPR:WrongInput','user name must be in the format name.surname.');
end

if contains(user, ' ')
    error('FSDA:ConfigClusterUNIPR:WrongInput','user name must be in the format name.surname without spaces.');
end


MATLABrelease = ['R' version('-release')];

def = struct;
def.MATLABrelease=MATLABrelease;

MapNetworkDrive='H:';

options = struct('NumWorkers',4, 'NumThreads', 2, ...
    'tasksxnode', 1, 'memxcpu','2G', 'IPaddress', true,...
    'MapNetworkDrive', MapNetworkDrive);

% tasksxnode may be absent for default option



if ~isempty(varargin)

    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:ConfigClusterUNIPR:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)

        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end


    def.user = user;
    def.NumWorkers = num2str(options.NumWorkers);
    def.NumThreads = options.NumThreads;
    def.tasksxnode = options.tasksxnode;
    def.memxcpu = options.memxcpu;    % string!
    def.MapNetworkDrive=options.MapNetworkDrive;

    if options.IPaddress ~= false
        if options.IPaddress == true
            tmp=char(java.net.InetAddress.getLocalHost.toString);
            tmp2=strsplit(tmp,'/');
            def.localIP= tmp2{2};
        else
            % TODO: add code for IP format validation
            def.localIP=options.IPaddress;
        end

    else
        % IPaddress == false, do not run UNIPRruncluster
        def.localIP=NaN;
    end
else
    % default arguments
    def.NumWorkers = num2str(4);
    def.NumThreads   = 2;
    def.tasksxnode = 1;
    def.memxcpu = '2G';
    tmp=char(java.net.InetAddress.getLocalHost.toString);
    tmp2=strsplit(tmp,'/');
    def.localIP= tmp2{2};
    def.MapNetworkDrive=MapNetworkDrive;
end

def.LocalJobStorageLocation= '';
def.JobStorageLocationOnPC = '';


% HPC unipr cluster parameters
def.user = user;
def.Type = 'remote';
def.ClusterMatlabRoot= [MATLABrelease ':/hpc/share/applications/matlab/' MATLABrelease];
def.ClusterHost= 'login.hpc.unipr.it';
def.RemoteJobStorageLocation = '/hpc/home';



profile = ['cluster ' MATLABrelease];

% Delete the old profile (if it exists)
ps = parallel.Settings;
pnidx = strcmp({ps.Profiles.Name},profile);
ws = warning('off');
ps.Profiles(pnidx).delete
warning(ws)


[~, hostname] = system('hostname');
hostname = strtrim(hostname);


tmp=strsplit(def.ClusterMatlabRoot,':');
def.ClusterMatlabRoot = tmp{2};


rootd = [def.LocalJobStorageLocation def.user];
loc = '.matlab';


jsl = fullfile(rootd,loc,'3p_cluster_jobs',hostname,MATLABrelease,def.Type);

if exist(jsl,'dir')==false
    [status,err,eid] = mkdir(jsl);
    if status==false
        error(eid,err)
    end
end


rootd = [def.RemoteJobStorageLocation def.user];


rjsl = [rootd '/' '.matlab' '/' '3p_cluster_jobs' '/' hostname '/' MATLABrelease '/' def.Type];

assembleClusterProfile(jsl, rjsl, def.user, profile, def);

pause(2);

if  ~isnan(def.localIP)
    UNIPRrunCluster(def.localIP);
end

end



function assembleClusterProfile(jsl, rjsl, user, profile, def)

% Create generic cluster profile
c = parallel.cluster.Generic;

c.IntegrationScriptsLocation = ['C:\ProgramData\MATLAB\SupportPackages\' def.MATLABrelease '\parallel\slurm\remote'];
c.NumWorkers = str2num(def.NumWorkers); %#ok<ST2NM>
c.OperatingSystem = 'unix';

% Depending on the submission type, populate cluster profile fields
% Set common properties for nonshared and remote
c.AdditionalProperties.Username = user;
c.AdditionalProperties.ClusterHost = def.ClusterHost;
c.ClusterMatlabRoot = def.ClusterMatlabRoot;

jsl = struct('windows',jsl,'unix',rjsl);

c.HasSharedFilesystem = true;
c.JobStorageLocation = jsl;


% UNIPR options before saving
clustuserpath= ['/hpc/home/' user '/parallel/matlab'];


localuserpath=[def.MapNetworkDrive '\parallel\matlab'];

c.JobStorageLocation = struct('windows', localuserpath, 'unix', clustuserpath);
% c.JobStorageLocation = struct('windows','H:\parallel\matlab', 'unix', clustuserpath);


winuser= getenv('username');
c.AdditionalProperties.IdentityFile=['C:\users\' winuser '\.ssh\id_rsa'];
c.AdditionalProperties.UseIdentityFile=true;
c.AdditionalProperties.IdentityFileHasPassphrase=false;
c.HasSharedFilesystem = true;
c.NumThreads= def.NumThreads;

% --partition=cpu,knl --tasks-per-node=1 --mem-per-cpu=2G
c.AdditionalProperties.AdditionalSubmitArgs = ...
    ['--partition=cpu,knl --tasks-per-node=' num2str(def.tasksxnode) ' --mem-per-cpu=' ...
    def.memxcpu];


% Save Profile
c.saveAsProfile(profile);
c.saveProfile('Description', profile)

% Set as default profile
parallel.defaultClusterProfile(profile);

end

