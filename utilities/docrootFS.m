function docr=docrootFS()
%docrootFS calls docroot or subfolder of FSDA according to Matlab version

a=ver('matlab');

if str2double(a.Version)>7.14
    docr=docroot;
else
        fsep=filesep;

    try
        FileWithFullPath=which('docsearchFS.m');
        [docr]=fileparts(FileWithFullPath);
        docr=[docr fsep 'helpfiles'];
    catch
        IndexFileFullPAth=[ docr fsep 'FSDA' fsep 'index.html'];
        web(IndexFileFullPAth);
    end
        docr=strrep(docr, '\', '/');

end

end
