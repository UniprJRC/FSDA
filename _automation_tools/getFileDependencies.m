function list =  getFileDependencies(file)
% searches the paths of the .m-Files the given file is using. This function
% depends on the depfun-Function of Matlab. Though it is not guaranteed,
% that all files will be found. E.g. files in evaluate-constructs won't be
% found!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters:
%     - file: the file for which the dependencies should be searched. This
%             could be only a string with the function-name, the filename
%             or the path to a file.
%
% Returns:
%     - list: a cell-Array with the paths to all files the given file uses
%             as strings. list is -1 if an error occured while processing.
%             See errormessage for details.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% author:   Torsten Hopp
%
% Copyright 2008-2019.
%
%$LastChangedDate::                      $: Date of the last commit

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% checking the input-parameter
if(isempty(file))
    disp('ERROR. Input Parameter is empty')
    list = -1;
    return;
end
if(~ischar(file))
    disp('ERROR. Input Parameter is not char')
    list = -1;
    return;
end


% list: a cell-array for the result
list = {which(file)};
% fileList: array of strings to search in the current iteration
fileList = list;

% While new files where found in last iteration
while (~isempty(fileList))
    fileListSize = length(fileList);
    newOnes = [];
    newCounter = 1;
    currentFileList = fileList;
    fileList = [];

    % running trough the current file list
    for k=1:fileListSize
        % calling depfun. -toponly indicates non-recursive search. -quiet
        % suppresses output to command line
        templist = depfun(currentFileList{k},'-toponly','-quiet');

        % the first file is already included 
        templist(1) = []; 
        
        % calling the helper-function to delete files which contain
        % <matlabroot> in the path from the list.
        templist = deleteMatlabFiles(templist);

        listSize = length(list);
        tempListSize = length(templist);

        % excluding first iteration with this if-clause
        if(listSize > 1)

            % comparing lists --> Searching for files that are not in the
            % list yet
            for i=1:tempListSize
                contained = 0;
                j=1;
                % first condition stops comparing if something was found.
                while contained == 0 & j<=listSize
                    if((strcmp(templist{i},list{j})==1) )
                        contained = 1;
                    end
                    j = j+1;
                end

                % if file isn't in the list yet...
                if(contained == 0)

                    % add it to a cellArray for new filenames...
                    % if-clause excludes first run with empty
                    % newOnes-variable
                    if(~isempty(newOnes))

                        % if not first run: comparing the filename with the
                        % filenames that are in the newOnes-List to avoid
                        % double-enties
                        a = strfind(newOnes,templist{i});
                        counter = 0;
                        for p=1:length(a)
                            if(~isempty(a{p}))
                                counter = counter +1;
                            end
                        end
                        % if not found in list: add it to the list!
                        if(counter==0)
                            newOnes{newCounter} = templist{i};
                            newCounter = newCounter + 1;
                        end
                    else
                        % if list is empty: add it without checking!
                        newOnes{newCounter} = templist{i};
                        newCounter = newCounter + 1;
                    end
                end
            end

            % if new filenames where found in this run...
            if(~isempty(newOnes))
                % ... search the new found files in the next iteration
                fileList = [fileList newOnes];
                % ... add the new found files to the list
                list(length(list)+1:length(list)+length(newOnes)) = newOnes;
                newOnes = [];
                newCounter = 1;
            end
        else
            % if first iteration: set the first-level dependecies as
            % filelist for next iteration and add them to the list.
            fileList = templist;
            list = [list; templist];
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function rList = deleteMatlabFiles(pList)
% deletes filenames with matlabroot in path from the given list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
%     - pList: a cell-Array of filenames as string
%
% Returns:
%     - rList: a cell-Array of filenames as string
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excludedDir = matlabroot;
indexes = strfind(pList,excludedDir);
keepIndexes = cellfun(@isempty,indexes);
rList = pList(keepIndexes);
