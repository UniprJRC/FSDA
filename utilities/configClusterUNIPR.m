function configClusterUNIPR(user, varargin)
% configClusterUNIPR create a cluster profile for MATLAB Parallel Server on
% the UNIPR HPC facilities and make the client Windows PC a member of the
% cluster.
%
% Using MATLAB an interactive parallel pool object on a Windows PC client
% connected to UNIPR HPC facilities is a complex task. To be able to
% exploit this feature the user need to create an ad hoc MATLAB Parallel
% Server profile, a valid ssh key pairs and, through SSH tunneling, insert
% the client Windows PC as a member of the HPC cluster.
%
%
%<a href="matlab: docsearchFS('configClusterUNIPR')">Link to the help function</a>
% Required input arguments:
%
%    user:      user SSH credentials of the HPC cluster. Character vector.
%               username specified as a char vector, in the format of name.surname
%               same as the mail account (without the domain part).
%
% Optional input arguments:
%
%
%  NumWorkers :  Indicator for constant term. Scalar. If 1, a model with
%               constant term will be fitted (default), else no constant
%               term will be included.
%               Example - 'intercept',1
%               Data Types - double
%
% NumThreads:   number of threads per CPU. Scalar. The default value is 2.
%               Example - 'NumThreads', 2
%               Data Types - double
%
%
%   tasksxnode: number of tasks per node. Scalar. The default value is 2.
%               Example - 'tasksxnode', 2
%               Data Types - double
%
%   memxcpu:    amount of RAM per CPU. String. The default value is '2G'.
%               Example - 'memxcpu','2G'
%               Data Types - double
%
% runclusterUNIPR:  initialize and connect the Windows client to the HPC cluster. Logical or string.
%               The default value is true.
%               Example - 'runclusterUNIPR', true
%               Data Types - boolean
%
% See also:     runclusterUNIPR
%
% Copyright 2008-2022.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('configClusterUNIPR')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%

%{
% all default options
configClusterUNIPR('john.doe')
%}

%{
% create or update the MATLAB parallel profile only
configClusterUNIPR('john.doe', 'runclusterUNIPR', false, 'tasksxnode', 1, 
'memxcpu','2G', 'NumWorkers', 2, 'NumThreads', 2)
%}

%{
% run also the function runclusterUNIPR() with client IP address
configClusterUNIPR('john.doe', 'runclusterUNIPR', true)
%}

%{
% run also the function runclusterUNIPR() with specific IP address
configClusterUNIPR('john.doe', 'runclusterUNIPR', '192.168.1.10')
%}


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

options = struct('NumWorkers',4, 'NumThreads', 2, ...
    'tasksxnode', 1, 'memxcpu','2G', 'runclusterUNIPR', true);

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

    if options.runclusterUNIPR ~= false
        if options.runclusterUNIPR == true
            tmp=char(java.net.InetAddress.getLocalHost.toString);
            tmp2=strsplit(tmp,'/');
            def.localIP= tmp2{2};
        else
            % TODO: add code for IP format validation
            def.localIP=options.runclusterUNIPR;
        end

    else
        % runclusterUNIPR == false, do not run runclusterUNIPR
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
end

% HPC unipr cluster parameters
def.user = user;
def.Type = 'remote';
def.ClusterMatlabRoot= [MATLABrelease ':/hpc/share/applications/matlab/' MATLABrelease];
def.ClusterHost= 'login.hpc.unipr.it';
def.LocalJobStorageLocation= '';
def.RemoteJobStorageLocation = '/hpc/home';
def.JobStorageLocationOnPC = '';


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
    runClusterUNIPR(def.localIP);
end

end



function assembleClusterProfile(jsl, rjsl, user, profile, def)

% Create generic cluster profile
c = parallel.cluster.Generic;


c.IntegrationScriptsLocation = 'C:\ProgramData\MATLAB\SupportPackages\R2021b\parallel\slurm\remote';
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
c.JobStorageLocation = struct('windows', 'H:\parallel\matlab', 'unix', clustuserpath);

winuser= getenv('username');
c.AdditionalProperties.IdentityFile=['C:\users\' winuser '\.ssh\id_rsa'];
c.AdditionalProperties.UseIdentityFile=true;
c.AdditionalProperties.IdentityFileHasPassphrase=false;
c.HasSharedFilesystem = true;
c.NumThreads= def.NumThreads;

% --partition=cpu,knl --tasks-per-node=1 --mem-per-cpu=2G
c.AdditionalProperties.SubmitArguments = ...
    ['--partition=cpu,knl --tasks-per-node=' num2str(def.tasksxnode) ' --mem-per-cpu=' ...
    def.memxcpu];


% Save Profile
c.saveAsProfile(profile);
c.saveProfile('Description', profile)

% Set as default profile
parallel.defaultClusterProfile(profile);

end

