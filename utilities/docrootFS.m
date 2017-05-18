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
        FSDAroot=FullPath(1:end-length(FileName)-3);
        FSDAhelproot=[FSDAroot filesep '_tmp_helpfiles' filesep 'FSDA'];
        
        sprintf('----------------------------------------------------------')
        fprintf('%s\n\n',['FSDA HTML help files must be copied inside path "' destinationpathFSDAdoc '"']);
        msg1=['The FSDA help files have not been copied under the MATLAB help folder: "' docroot '".'];
        msg2='  Please check the user right permissions and possibly move manually the folder "';
        msg3='"  in the MATLAB help folder. For additional details please see file: ';
        msg4=[FSDAroot filesep 'InstallationNotes.pdf'];
        
        msg=regexprep([msg1 msg2 FSDAhelproot msg3 msg4],'\\','/');
        warning('FSDA:docrootFS:wrongSetUp',msg)
        h=mydialog(msg);
        import com.mathworks.mlwidgets.html.HTMLRenderer;
        % create component
        r = HTMLRenderer;
        % set the text to display
        r.setHtmlText(['<html> <a href="' msg4 '">Installation Notes</a></html>']);
        % make sure the component is opaque
        r.setOpaque(true);
        % add the component
        javacomponent(r, [120 49 100 30], h);
        %[left, bottom, width, height]
        
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

function d=mydialog(msg)
d = dialog('Position',[300 300 380 250],'Name','FSDA warning message');

uicontrol('Parent',d,...
    'Style','text',...
    'Position',[20 80 330 155],...
    'String',msg,'FontSize',11);

uicontrol('Parent',d,...
    'Position',[135 20 70 25],...
    'String','Close',...
    'Callback','delete(gcf)');
end