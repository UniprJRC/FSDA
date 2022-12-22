%% Script that automatically adds description to the newest datasets of FSDA

clear

load alldescriptions.mat %mat file that contains all the descriptions

% Get a list of all variables in the workspace
varnames = whos;


% Loop through the variables and check if they are structs (only the
% descriptions are structs in here so basically it is extracting every
% description in the mat file, without considering other variables that may
% have slipped in by mistake when creating the .mat file with the
% descriptions

n = numel(varnames);
names=cell(n,1);
for j= 1:n
    % verify if the variable is a 'struct'
    if strcmp(varnames(j).class, 'struct')
        names{j}=varnames(j).name;
    end
end

%% this part of the code is tasked to store the descriptions in a struct for
% easier indexing in the next steps
descr=struct;
for j= 1:n
    if ~isempty(names{j}) %verify that the 'i' element in the array is not empty
        descr.(names{j}) = eval(names{j});
    end
end


%% 'folder' is the current folder
Folder= pwd;

%dir looks in the directory contained in 'Folder', while fullfile creates
% a full file name from the separate objects i give to it
Filelist= dir(fullfile(Folder,'**/*.mat')); %creating a list of all the datasets to edit
Nfile= numel(Filelist);
pat=".mat";


%% for loop that should set descriptions for each file

for j=1:Nfile
    Filej= Filelist(j).name;
    name= extractBefore(string(Filej), pat);
    %check if the file to be loaded has a description stored, so that the
    %script continues running. Warning: it may leave behind some files that
    %are missing a description, if their description is misspelled!
    if isfield(descr, name)
        Filepathj= Filelist(j).folder;
        dirFull= fullfile(Filepathj, Filej);
        tmp= load(dirFull);
        tmp.(name).Properties.Description = descr.(name).Properties.Description;
        save(Filej, '-struct', 'tmp', name) %saves the new file in current folder
    end
end

%% Show the file .mat which still do not have a description
clear
Folder=pwd;
Filelist= dir(fullfile(Folder,'**/*.mat')); %creating a list of all the datasets to edit
Nfile= numel(Filelist);
pat=".mat";

for j=1:Nfile
    Filej= Filelist(j).name;
    name= extractBefore(string(Filej), pat);
    if isempty(intersect(name,{'alldescriptions' 'facemasksALL'}))
        eval(['tmp=load(''' char(name) ''');']);
        if  isempty(tmp.(name).Properties.Description)
            disp(Filej)
        end
    end
end

