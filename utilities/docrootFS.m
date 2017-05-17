function docr=docrootFS()
%docrootFS calls docroot or subfolder of FSDA according to Matlab version

a=ver('matlab');

if str2double(a.Version)>7.14
    docr=docroot;
    destinationpathFSDAdoc=[docr filesep 'FSDA'];
    chkFSDAinsideMatlab=exist(destinationpathFSDAdoc,'file');
    if chkFSDAinsideMatlab ~=7
        FileName='addFSDA2path';
        FullPath=which(FileName);
        FSDAhelproot=[FullPath(1:end-length(FileName)-3) filesep 'helpfiles' filesep 'FSDA'];
        
        sprintf('----------------------------------------------------------')
        fprintf('%s\n\n',['FSDA HTML help files must be copied inside path ' destinationpathFSDAdoc]);
        msg1=['The FSDA help files have not been copied under the MATLAB help folder: ' docroot '.'];
        msg2='  Please check the user right permissions and possibly copy manually the folder  ';
        msg3='  in the MATLAB help folder';
        msg=regexprep([msg1 msg2 FSDAhelproot msg3],'\\','/');
        warning('FSDA:docrootFS:wrongSetUp',msg)
        mydialog(msg)
        error('FSDA:docrootFS:wrongSetUp','Please read warning above and take appropriate action!')
    else
    end
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

function mydialog(msg)
    d = dialog('Position',[300 300 350 150],'Name','FSDA warning message');

     uicontrol('Parent',d,...
               'Style','text',...
               'Position',[20 80 310 60],...
               'String',msg);

     uicontrol('Parent',d,...
               'Position',[135 20 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');
end