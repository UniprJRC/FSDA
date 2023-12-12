function out=gitKeysSSHroot
% addGitSshRsaKeys adds the user's public RSA keys to the
% local git repository to avoid continuos request of GitHub personal token
% at each push to Origin

if ispc
    username=getenv('USERPROFILE');
else
    username=getenv('HOME');
end
try
    git = settings().matlab.sourcecontrol.git;
    git.PrivateKeyFile.PersonalValue = [username '/.ssh/id_rsa'];
    git.PublicKeyFile.PersonalValue = [username '/.ssh/id_rsa.pub'];
catch ME
    disp('error',ME)
end

out = struct;
out.PrivateKeyFile=git.PrivateKeyFile.PersonalValue;
out.PublicKeyFile=git.PublicKeyFile.PersonalValue;
if nargout == 0
    disp(out)
end

end