function UNIPRrunCluster(IPaddress)
%UNIPRrunCluster  connects local machine to a remote cluster for parallel computing
%
%<a href="matlab: docsearchFS('UNIPRrunCluster')">Link to the help function</a>
%
%  This routine assumes that at least a cluster profile on a remote server has
%  already been prepared and enables to connect the local computer to the
%  cluster. The cluster profile (together with the settings of the SLURM
%  parameters) can be created with function UNIPRconfigCluster.m. 
%  Once UNIPRrunCluster is launched a pop up window automatically appears
%  which enables the user to select the cluster to which local computer has
%  to connect. This function can be called without input arguments. In this
%  case IPaddress of the local machine is used. Alternatively the first
%  input argument specifies the IP address which has to be used.
%
%
% Required input arguments:
%
%
% Optional input arguments:
%
%    IPaddress   : IPV4 address of the machine which has to be connected to
%                  the remote cluster. Character or string.
%                  A valid IPV4 address has a format xxx.xxx.xxx.xxx where
%                  xxx has a range from 0 to 255. Note that the
%                  IPaddress which is supplied is not validated.
%               Example - '192.168.178.22' or "160.78.248.185"
%               Data Types - double
%
%
% Output:
%
%
% See also: UNIPRconfigCluster, pctconfig
%
% References:
%
%
% Acknowledgements:
% 
% This function is a modification of a function which had
% originally been written by Fabio Spataro <fabio.spataro@unipr.it> of the
% High Performance Computing Facilities of the University of Parma.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('UNIPRrunCluster')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{ 
    Call to UNIPRrunCluster without IP address.
    UNIPRrunCluster
%}

%{
    Call to UNIPRrunCluster with IP address in character format.
    UNIPRrunCluster('160.78.6.11')
%}

%{
    Call to UNIPRrunCluster with IP address in string format.
    UNIPRrunCluster("160.78.6.11")
%}


%% Beginning of code

% MATLAB version check
if verLessThanFS([9 11])
    error('FSDA:ConfigClusterUNIPR:WrongMATLABVersion','At least MATLAB R2021b is needed.');
end


if nargin == 0  % use current IP address
    tmp=char(java.net.InetAddress.getLocalHost.toString);
    tmp2=strsplit(tmp,'/');
    IPaddress = tmp2{2};
else
    IPaddress=char(IPaddress);
end

% Cluster profiles
[list, default] = parallel.clusterProfiles;

% Default profile index
default_indx = find(string(list) == default);

% List selection dialog box
[indx, tf] = listdlg(...
    'PromptString', 'Select a cluster profile:',...
    'SelectionMode', 'single',...
    'ListSize', [300 100],...
    'InitialValue', default_indx,...
    'ListString', list);

% Return if "Cancel button" is selected
if tf == 0
    return
end

% Cluster profile
cluster = parcluster(string(list(indx)));

% Gets properties from the cluster profile
try
    Username     = cluster.AdditionalProperties.Username;
    ClusterHost  = cluster.AdditionalProperties.ClusterHost;
    IdentityFile = cluster.AdditionalProperties.IdentityFile;
catch
    warning('Problem using profile ''' + string(list(indx)) + '''.');
    return
end

% Port range
rng shuffle;
port1 = randi([49152, 65535 - 1], 1, 1);
port2 = port1 + 1;
PortRange = [port1 port2];

% Command formats
tunnel_cmd_format = "ssh %s@%s -i %s -N -R %s:%d:%s:%d";
monitor_cmd_format = "ssh %s@%s -i %s TERM=xterm-color watch sacct --allocations --user=%s";

% Commands
cmd1 = sprintf(tunnel_cmd_format, Username, ClusterHost, IdentityFile, ClusterHost, port1, IPaddress, port1);
cmd2 = sprintf(tunnel_cmd_format, Username, ClusterHost, IdentityFile, ClusterHost, port2, IPaddress, port2);
cmd3 = sprintf(monitor_cmd_format, Username, ClusterHost, IdentityFile, Username);

% Operating system dependent commands
if ismac
    % Code to run on Mac platform
    oscmd1 = sprintf("xterm -geometry 120x10 -bg 'white' -title '%s' -e '%s' &", cmd1, cmd1);
    oscmd2 = sprintf("xterm -geometry 120x10 -bg 'white' -title '%s' -e '%s' &", cmd2, cmd2);
    oscmd3 = sprintf("xterm -geometry 120x40 -bg 'white' -title '%s' -e '%s' &", cmd3, cmd3);
elseif isunix
    % Code to run on Linux platform
    oscmd1 = sprintf("xterm -geometry 120x10 -bg 'white' -title '%s' -e '%s' &", cmd1, cmd1);
    oscmd2 = sprintf("xterm -geometry 120x10 -bg 'white' -title '%s' -e '%s' &", cmd2, cmd2);
    oscmd3 = sprintf("xterm -geometry 120x40 -bg 'white' -title '%s' -e '%s' &", cmd3, cmd3);
elseif ispc
    % Code to run on Windows platform
    oscmd1 = sprintf('start "%s" %s', cmd1, cmd1);
    oscmd2 = sprintf('start "%s" %s', cmd2, cmd2);
else
    disp('Platform not supported')
    return
end

% Executes operating system dependent commands
[status1,results1] = system(oscmd1);
[status2,results2] = system(oscmd2);
if (ismac || isunix)
    [status3,results3] = system(oscmd3);
end

% Removes some variables from the currently active workspace
clearvars indx tf
clearvars client_addr
clearvars port1 port2
clearvars tunnel_cmd_format
clearvars monitor_cmd_format
clearvars cmd1
clearvars cmd2
clearvars cmd3
clearvars oscmd1
clearvars oscmd2
if (ismac || isunix)
    clearvars oscmd3
end

% Configure a Parallel Computing Toolbox session
pctconfig('hostname', ClusterHost);
pctconfig('portrange', PortRange);
end
