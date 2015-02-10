function Create_automatic_links()
%Create_automatic_links puts links in html documentation files 
%
% The purpose of this file is to automatically create the links to
% automatically navigate in alfabethical order among the html pages that
% are create. 
% REMARK: this routine scans current directory to find HTML files and put
% all alphabetical links in all the pages except for the files named
% users_guide.html
% function_reference.html
% alphabetical.html
% acknowledgments.html
% bibliography.html
% fsda_product_page.html
% fsda_toolbox_product_page.html
% datasets.html
% exampleindex.html
% developers.html
% release_notes.html
%  datasets_clu
%  datasets_mv
%  datasets_reg
%  cluster_intro
%
% Similarly all .html files which start with the following words
% mult
% transf
% regression
% statistical
% group
% automatic_outl_detect_fsm
%
% are ignored

% Copyright 2008-2015.
% Written by FSDA team
%
% Last modified 06-Feb-2015

% Examples:

%
%{
    % Navigate in a folder which contains the html files and then digit
    Create_automatic_links()
%}

%% Beginning of code
% List all file with .html extension in current directory
namefiles=dir('*.html');


% Find how many they are
totnum=size([namefiles.isdir],2);

%% Put all of them into a cell array of strings
filenames=cell(totnum,1);
% filenamesno contains all the files for which links do not have to be
% created
filenamesno=cell(totnum,1);
ij=1;
ino=1;
for i=1:totnum
    if strncmp('intro',namefiles(i).name,5)+ strncmp('product',namefiles(i).name,7) ...
            + strncmp('users_guide',namefiles(i).name,11) ...
            + strncmp('function_reference',namefiles(i).name,18)  ...
            + strncmp('alphabetical',namefiles(i).name,12)  ...
            + strncmp('acknowledgments',namefiles(i).name,15)  ...
            + strncmp('bibliography',namefiles(i).name,12)  ...
            + strncmp('mult',namefiles(i).name,4)  ...
            + strncmp('transf',namefiles(i).name,6)  ...
            + strncmp('regression',namefiles(i).name,10)  ...
            + strncmp('fsda_product_page',namefiles(i).name,17)  ...
            + strncmp('fsda_toolbox_product_page',namefiles(i).name,25)  ...
            + strncmp('statistical',namefiles(i).name,11)  ...
            + strncmp('exampleindex',namefiles(i).name,12)  ...
            + strncmp('developers',namefiles(i).name,10) ...
            + strncmp('group',namefiles(i).name,5) ...
            + strncmp('datasets',namefiles(i).name,8) ...
            + strncmp('datasets_clu',namefiles(i).name,12) ...
            + strncmp('datasets_mv',namefiles(i).name,11) ...
            + strncmp('datasets_reg',namefiles(i).name,12) ...
            + strncmp('cluster_intro',namefiles(i).name,13) ...
            + strncmp('automatic_outl_detect_fsm',namefiles(i).name,25) ...
            + strncmp('release_notes',namefiles(i).name,13) ==0
        
        
        filenames{ij} = namefiles(i).name;
        ij=ij+1;
    else
        filenamesno{ino} = namefiles(i).name;
        ino=ino+1;
    end
end
filenames(ij-1)
% sort them (from A to z)
filenames=sort(filenames(1:ij-1));
filenamesno=filenamesno(1:ino-1);
disp(filenamesno)

%% LOOP through all filenames (excluding filenamesno) to create all the
% automatic links

% ii=7;

for ii= 2:(size(filenames,1)-1)
    % Access and open the ith file
    filen=filenames(ii);
    fileID = fopen(char(filen), 'r+');
    
    % fileID = fopen(filen, 'r');
    % all = textscan(fileID, '%s');
    
    % Insert the file into fstring
    fstring=fscanf(fileID,'%c');
    
    % Replace right link on top of the page
    
    % Find the arrow RIGHT link on TOP of the page
    POSbprev=regexp(fstring,'b_prev.gif');
    POSinilink=regexp(fstring,'<a href=');
    POSinilink=POSinilink(POSinilink>POSbprev(1));
    POSinilink=POSinilink(1)+9;
    
    POSendlink=regexp(fstring,'">');
    POSendlink=POSendlink(POSendlink>POSinilink)-1;
    
    % string which has to be replaced
    % fstring(POSinilink:POSendlink)
    
    % Replace the old link with the new link
    fstring=[fstring(1:POSinilink-1)  char(filenames(ii+1)) fstring(POSendlink+1:end)];
    
    % Write fstring into a file
    % count=fprintf(fileID,'%s',fstring);
    
    % Replace left link on top of the page
    
    % Find the arrow LEFT link on TOP of the page
    POSinilink=regexp(fstring,'<a href=');
    POSinilink=POSinilink(1)+9;
    
    POSendlink=regexp(fstring,'">');
    POSendlink=POSendlink(POSendlink>POSinilink)-1;
    
    % arrow LEFT link on TOP of the page which has to be replaced
    % fstring(POSinilink:POSendlink)
    
    % Replace the old link with the new link
    fstring=[fstring(1:POSinilink-1)  char(filenames(ii-1)) fstring(POSendlink+1:end)];
    
    
    % Replace links and labels at the bottom of the page
    
    % link on the LEFT arrow at the BOTTOM of the page
    % Find the second instance of b_prev.gif
    POSbprev=regexp(fstring,'b_prev.gif');
    POSbprev=POSbprev(2);
    
    POSinilink=regexp(fstring,'<a href=');
    POSinilink=POSinilink(POSinilink<POSbprev);
    POSinilink=POSinilink(end)+9;
    
    POSendlink=regexp(fstring,'">');
    POSendlink=POSendlink(POSendlink>POSinilink)-1;
    
    % string which has to be replaced
    % fstring(POSinilink:POSendlink)
    
    % Replace the old link with the new link
    fstring=[fstring(1:POSinilink-1)  char(filenames(ii-1)) fstring(POSendlink+1:end)];
    
    % Replace the label close to the LEFT arrow at the BOTTOM of the page
    % Find </td>
    POSendlink=regexp(fstring,'</td>');
    % Among all instances greater than bprev take the second
    POSendlink=POSendlink(POSendlink>POSinilink);
    POSendlink=POSendlink(2)-1;
    
    % Find string >
    % the position is the greatest among those smaller than POSendlink
    POSinilink=regexp(fstring,'>');
    POSinilink=POSinilink(POSinilink<POSendlink);
    POSinilink=POSinilink(end)+1;
    % string which has to be replaced
    % fstring(POSinilink:POSendlink)
    
    % Replace the old string with the new string
    fstring=[fstring(1:POSinilink-1)  char(filenames(ii-1))  fstring(POSendlink+1:end)];
    
    
    
    % link on the RIGHT arrow at the BOTTOM of the page
    
    % Find the second instance of b_next.gif
    POSbnext=regexp(fstring,'b_next.gif');
    POSbnext=POSbnext(2);
    
    POSinilink=regexp(fstring,'<a href=');
    POSinilink=POSinilink(POSinilink<POSbnext);
    POSinilink=POSinilink(end)+9;
    
    POSendlink=regexp(fstring,'">');
    POSendlink=POSendlink(POSendlink>POSinilink)-1;
    
    % string which has to be replaced
    % fstring(POSinilink:POSendlink)
    
    % Replace the old link with the new link
    fstring=[fstring(1:POSinilink-1)  char(filenames(ii+1))  fstring(POSendlink+1:end)];
    
    % Replace the label close to the RIGHT arrow at the BOTTOM of the page
    % Find </td>
    POSendlink=regexp(fstring,'</td>');
    % Among all instances smaller than bprev take the second largest
    POSendlink=POSendlink(POSendlink<POSinilink);
    POSendlink=POSendlink(end)-1;
    
    % Find string >
    % the position is the greatest among those smaller than POSendlink
    POSinilink=regexp(fstring,'>');
    POSinilink=POSinilink(POSinilink<POSendlink);
    POSinilink=POSinilink(end)+1;
    % string which has to be replaced
    % fstring(POSinilink:POSendlink)
    
    % Replace the old string with the new string
    fstring=[fstring(1:POSinilink-1)  char(filenames(ii+1)) fstring(POSendlink+1:end)];
    
    fclose('all');
    
    % file1ID=fopen('newfile.html','w');
    % Reopen the file for writing
    file1ID=fopen(char(filen),'w');
    fprintf(file1ID,'%s',fstring);
    fclose('all');
    
end
end
