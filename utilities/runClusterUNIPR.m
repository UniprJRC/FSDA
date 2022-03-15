function runClusterUNIPR(client_addr)
%PCTINIT Initialize Parallel Computing Toolbox.
%	pctinit(client_addr)
%
%	EXAMPLES:
%	pctinit('192.168.178.22');
%	pctinit('192.168.178.24');
%	pctinit('160.78.248.185');
%
%	INPUT:
%	Client_addr - Client IP address.
%
%	Fabio Spataro <fabio.spataro@unipr.it>, 2021-11-02

% MATLAB vesrion check
if verLessThanFS([9 11])
    error('FSDA:ConfigClusterUNIPR:WrongMATLABVersion','At least MATLAB R2021b is needed.');
end

% Number input arguments
if nargin == 0
    tmp=char(java.net.InetAddress.getLocalHost.toString);
    tmp2=strsplit(tmp,'/');
    client_addr = tmp2{2};
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

% Return if "Cancel" is selected
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
cmd1 = sprintf(tunnel_cmd_format, Username, ClusterHost, IdentityFile, ClusterHost, port1, client_addr, port1);
cmd2 = sprintf(tunnel_cmd_format, Username, ClusterHost, IdentityFile, ClusterHost, port2, client_addr, port2);
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
    return;
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
