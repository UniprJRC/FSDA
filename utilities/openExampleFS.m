function openExampleFS(filename)
%openExampleFS Open an example for modification and execution.
    examp=which(filename);
    
    examp1=strrep(examp,'\','\\');
    opentoline(examp1,5)
end
