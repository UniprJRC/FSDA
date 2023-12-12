function c=setdiffFS(a,b)
%setdiffFS finds the positive integers in a which are not present in the positive integers in b
%
% This is a faster special case of function setdiff when
% both vectors in a and b just contain positive integer numbers.
%
%<a href="matlab: docsearchFS('setdiffFS')">Link to the help function</a>
%
% Required input arguments:
%
%    a:         vector containing positive integer elements. Vector. A
%               vector of length na containing positive integer numbers.
%
%    b:         vector containing positive integer elements. Vector. A
%               vector of length nb containing positive integer numbers.
%
%
% Optional input arguments:
%
%
% Output:
%
%    c:         vector containing positive integer elements thare are on a but not in b.
%               Column vector. 
%               Note that the elements of c contain no repetitions and are sorted.
%
%
% See also: setdiff
%
% References:
%
%
% Riani, M., Perrotta, D. and Cerioli, A. (2015), The Forward Search for
% Very Large Datasets, "Journal of Statistical Software"
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('setdiffFS')">Link to the help page for this function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of use of setdiffFS.
    % Define two vectors (containing positive integers) with values in
    % common.
    A = [3 6 2 1 5 1 1]; 
    B = [2 4 6];
    C=setdiffFS(A,B);
    disp(C);
%}

%{
    % Time comparison with setdiff.
    % 20000 calls to setdiff and to setdiffFS.
    % Analysis of computational time.
    n=100;
    nsimul=20000;
    tSETDIFF=0;
    tSETDIFFFS=0;
    for j=1:nsimul
        a=randi(n,[300,1]);
        b=randi(n,[40,1]);

        tsetdiff = tic;
        c=setdiff(a,b);
        tSETDIFF = tSETDIFF + toc(tsetdiff);

        tsetdiffFS = tic;
        cFS=setdiffFS(a,b);
        tSETDIFFFS = tSETDIFFFS + toc(tsetdiffFS);

        if ~isequal(c,cFS)
            error('FSDA:setdiffFS:WrongOutput','c and cFS are different')
        end

    end

    disp(array2table([tSETDIFF tSETDIFFFS],'VariableNames',{'setdiff time' 'setdiffFS time'}))
%}

%% Beginning of code

ma    = max([a(:);b(:)]);
aT    = false(ma,1);
bT    = aT;
aT(a) = true;
bT(b) = true;

% c = vector containing numbers which are inside vector a which are not
% present in b. Elements in c are sorted and contain no repetitions.
c = find(bitxor(aT , bitand(aT,bT)));
% The instruction above is faster than the one below
% c = find(bitand(aT , not(bT)));
% The instruction above is faster than the one below
% c = find(aT & ~bT);

end
%FScategory:UTIGEN