%%import descriptions from xlsx file and create an alldescriptions.mat file

clear

filename= 'Description1.xlsx'; %the file that contains the descriptions
data= readcell(filename);
%this part removes any missing lines in the descriptions as it breaks the
%code later on
dn= data(:,1);
dd= data(:,2);
rmis= ~cellfun(@(x) any(isa(x,'missing')), dd);
desc= cat(2, dn(rmis), dd(rmis));
d1=desc(:,1);

%name contains all the name of the description, removing empty lines
mask= ~cellfun(@(x) any(isa(x,'missing')), d1);
name= d1(mask);

%% 'while' loop that creates description variables
descr= struct();
l= numel(desc(:,2));
j=1;
while j<l
    namJ= string(desc(j,1));
    %the next if statement detects if there is only a single line of
    %description or multiple.
    if ~isempty(intersect(namJ, name))
        descr.(namJ).Properties.Description = string(desc(j,2));
        j= j+1;
        %if there are multiple lines of description 'vertcat' concatenates
        %them vertically, as they are written in the FSDA documentation
    else
        n1= string(desc(j-1,1));
        i=j;
        t=i-1;
        D = vertcat(string(desc(t,2)), string(desc(i,2)));
        %this next if statement prevents the code from breaking if there
        %are only 2 lines of description. otherwise it would sometimes bug
        %and add the description of the following one
        if ~isempty(intersect(string(desc(i+1,1)), name))
            j= i+1;
            continue 
        else
            while mask(j) == 0
                i = i+1;
                D = vertcat(D, string(desc(i,2)));
                %the next if statement prevents the 'Index in position 1
                %exceeds array bounds.' Error from occurring, and breaks the
                %nested while loop
                if i>l-1
                    j=i+1;
                    break
                    %the next if statement detects if the following line is a
                    %new description or part of the old one
                elseif ~isempty(intersect(string(desc(i+1,1)), name))
                    j= i+1;
                    continue
                end
            end
        end
        %when all the strings composing the description have been
        %concatenated, the new description is saved in the struct
        D = strjoin(D, '\n');
        descr.(n1).Properties.Description = D;
    end
end

%% save the descriptions in the .mat file 'alldescriptions.mat'

save("alldescriptions.mat", '-struct', 'descr')