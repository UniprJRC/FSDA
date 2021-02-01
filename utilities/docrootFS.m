function docr=docrootFS()
%docrootFS calls docroot or subfolder of FSDA according to Matlab version

% Copyright 2008-2021.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit


a=ver('matlab');

if str2double(a.Version)>7.14
    docr=docroot    ;
    destinationpathFSDAdoc=[docr filesep 'FSDA'];
    chkFSDAinsideMatlab=exist(destinationpathFSDAdoc,'file');
    if chkFSDAinsideMatlab ~=7
        FileName='addFSDA2path';
        FullPath=which(FileName);
        FSDAroot=FullPath(1:end-length(FileName)-3);
        FSDAhelproot=[FSDAroot filesep 'helpfiles' filesep 'FSDA'];
        
        char10=char(10);
        sprintf('----------------------------------------------------------')
        fprintf('%s\n\n',['FSDA HTML help files which are inside path' char10  '"' FSDAhelproot '"']);
        fprintf('%s\n\n',['must be copied inside MATLAB doc path' char10  '"' destinationpathFSDAdoc '"']);
        
        msg1=['The FSDA help files have not been copied under the MATLAB doc folder:' char10];
        msg1bis=['"' docroot '".' char10];
        msg2=['Please check the user right permissions and possibly copy manually the folder' char10 ];
        msg2bis= ['"' FSDAhelproot '"' char10];
        msg3=['in the MATLAB help folder.' char10];
        msg3bis= ['"' destinationpathFSDAdoc '"' char10];
        msg3tris=['Alternatively, run routine' char10 char 10 'installHelpFiles.m.' char10 char10];
        msg3fourth=['or click on the Copy files button below' char10 char10 'For additional details please see section "FSDA documentation" of file:' char10];
        msg4=[FSDAroot filesep 'doc' filesep 'GettingStarted.mlx'];
        
        msg=regexprep([msg1 msg1bis msg2 msg2bis msg3 msg3bis msg3tris msg3fourth msg4],'\\','/');
        warning('FSDA:docrootFS:wrongSetUp',msg)
        %h=mydialog(msg);
        answer = questdlg(msg,'Copy FSDA HTML help files','Copy files','Close','Close');
        drawnow; pause(0.05);
        if isequal(answer,'Copy files')
            installHelpFiles();
            if ~exist('docr','var')
                docr=docroot;
            end
            return
        else
            docr='http://rosa.unipr.it/';
            return
        end

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
% 
% function d=mydialog(msg)
% % Dimensione overall dialogue box
% d = dialog('Resize','on','Position',[300 300 580 400],'Name','FSDA warning message');
% 
% uicontrol('Parent',d,...
%     'Style','text',...
%     'Position',[20 50 530 325],...
%     'String',msg,'FontSize',11);
% 
% % Dimension of the buttom close
% uicontrol('Parent',d,...
%     'Position',[135 20 70 25],...
%     'String','Close',...
%     'Callback','delete(gcf)');
% end
