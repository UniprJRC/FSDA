function out = rcontFS( I, J, nrowt, ncolt, varargin)
%rcontFS generates a random two-way table with given marginal totals
%
%
%<a href="matlab: docsearchFS('rcontFS')">Link to the help function</a>
%
%  Required input arguments:
%
%
%         I   :   Number of rows of the simulated contingency table.
%                 Scalar.
%                 Scalar which contains the requested number of rows the
%                 output matrix must have.
%         J   :   Number of columns of the simulated contingency table.
%                 Scalar.
%                 Scalar which contains the requested number of columns the
%                 output matrix must have.
%       nrowt  :  Row totals of the simulated contingency table.
%                 Vector.
%                 Vector of length I containing the requested row totals the
%                 output matrix must have. First element refers to the
%                 total number of elements in the first row, ..., $I$-th
%                 element refers to the total number of elements in the
%                 $I$-th row. In other words, $nrowt=(n_{1.}, n_{2.}, \ldots, n_{I.})$.
%       ncolt  :  Column totals of the simulated contingency table.
%                 Vector.
%                 Vector of length J containing the requested column totals the
%                 output matrix must have. First element refers to the
%                 total number of elements in the first column, ..., $J$-th
%                 element refers to the total number of elements in the
%                 $J$-th column.
%                 In other words, $ncolt=(n_{.1}, n_{.2}, \ldots, n_{.J})$.
%
%
%
%  Optional input arguments:
%
%    nocheck  : Checks on input arguments. Boolean.
%               If nocheck is false (default) program checks whether
%               1) nrow and ncol are greater than 1;
%               2) min(nrowt) and min(ncolt) are strictly greater than 0;
%               3) length(nrowt)=I, and  length(ncol)=J;
%               4) the sum of the elements of vector nrowt is equal to the sum of the
%                  elements of vector ncolt, (in other words, we check whether the row and
%                  column sum vectors have the same grand total).
%               To avoid all the above checks set nocheck to true.
%               Example - 'nocheck',true
%               Data Types - boolean
%
%  algorithm  : Algorithm to use to create the random contingency table.
%               Character.
%               Character which specifies which algorithm must be used to
%               create the contingency table.
%               Possible values for algorithm are:
%               '144' in  this case the algorithm due to Boyett (1979) is
%               used and the output structure will contain field
%               out.matrix144.
%               '159' in  this case the algorithm due to Patefield (1979) is
%               used and the output structure will contain field
%               out.matrix159.
%               'all' in  this case the algorithms due
%               to Boyett (1979) and to Patefield (1981) are
%               both used. The output structure out will contain both
%               out.matrix144 and out.matrix159.
%               Example - 'algorithm','144'
%               Data Types - character
%
%  Output:
%
%  out :     A structure containing the following fields
%
%
% 		out.m144      =   $I$-by-$J$-table containing contingency table
%                         generated using algorithm '144' due to Boyett.
%                         The $(i,j)$-th element is equal to $n_{ij}$,
%                         $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$.
%                         This field is present only if input option
%                         algorithm is 'all' or '144'.
%                         Note that sum(out.144,1)=nrowt and that
%                         sum(out.144,1)=ncolt.
%
% 		out.m159      =   $I$-by-$J$-table containing contingency table
%                         generated using algorithm '159' due to Patefield.
%                         The $(i,j)$-th element is equal to $n_{ij}$,
%                         $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$.
%                         This field is present only if input option
%                         algorithm is 'all' or '159'.
%                         Note that sum(out.159,1)=nrowt and that
%                         sum(out.159,1)=ncolt.
%
%
% See also: CorAna, crosstab
%
%
% References:
%
% Boyett, J. (1979), Algorithm AS 144: Random R x C Tables with Given Row
% and Column Totals, "Applied Statistics", Vol. 28, pp. 329-332.
% Patefield, M. (1981), Algorithm AS 159: An Efficient Method of Generating
% RXC Tables with Given Row and Column Totals, "Applied Statistics", 
% Vol. 30, pp. 91-97.
%
% Acknowledgements:
%
% This routine is based on the codes written in FORTRAN77 version by
% Michael Patefield and James Boyett and on the Matlab version by John
% Burkardt.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('rcontFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %% rcontFS with all the default options.
    % Generate a random contingency table with 2 rows and 3 columns.
    nrow=2;
    ncol=3;
    % Fix the marginals of the two rows
    nrowt=[20 30];
    % Fix the marginals of the three columns
    ncolt=[25 15 10];
    % Generate the contingency table
    out=rcontFS(nrow,ncol,nrowt,ncolt)
    % Random contingency table based on algorthm AS144
    disp(out.m144)
    % Random contingency table based on algorthm AS159
    disp(out.m159)
%}

%{
    % rcontFS with option nocheck set to true.
    % Generate a random contingency table with 2 rows and 3 columns and do
    % not check input arguments.
    nrow=2;
    ncol=3;
    % Fix the marginals of the two rows
    nrowt=[20 30];
    % Fix the marginals of the three columns
    ncolt=[25 15 10];
    % Generate the contingency table and avoid checks
    out=rcontFS(nrow,ncol,nrowt,ncolt,'nocheck',true)
    % Random contingency table based on algorthm AS144
    disp(out.m144)
    % Random contingency table based on algorthm AS159
    disp(out.m159)
%}


%% Beginning of code

algorithm='all';
nocheck=0;
options=struct('algorithm',algorithm,'nocheck',nocheck);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:rcontFS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    algorithm=options.algorithm;
    nocheck=options.nocheck;
end


%% Beginning of code
if nocheck == false
    if I <= 1
        error('FSDA:rcontFS:WrongInputOpt','Number of rows of the contingency table cannot be smaller than 1')
    end
    
    if J <= 1
        error('FSDA:rcontFS:WrongInputOpt','Number of columns of the contingency table cannot be smaller than 1')
    end
    
    if min(nrowt)<=0
        error('FSDA:rcontFS:WrongInputOpt','Row totals must be integers strictly greater than 0')
    end
    
    if min(ncolt)<=0
        error('FSDA:rcontFS:WrongInputOpt','Column totals must be integers strictly greater than 0')
    end
    
    if length(nrowt)~=I
        error('FSDA:rcontFS:WrongInputOpt',['length of vector which contains marginal row totals bust be equal to nrow.'...
            ' Now length(nrowt)=' num2str(length(nrowt)) ' and nrow='  num2str(I)])
    end
    
    if length(ncolt)~=J
        error('FSDA:rcontFS:WrongInputOpt',['length of vector which contains marginal row totals bust be equal to nrow.'...
            ' Now length(nrowt)=' num2str(length(ncolt)) ' and nrow='  num2str(J)])
    end
    
    if sum(ncolt)~=sum(nrowt)
        error('FSDA:rcontFS:WrongInputOpt',['column sum vectors have the same sum.'...
            ' Now sum(nrowt)=' num2str(sum(nrowt)) ' and sum(ncolt)='  num2str(ncoly)])
    end
end


% nsubt = partial cumulative sum of column counts
nsubt=cumsum(ncolt);
ntotal = nsubt(J);

out=struct;

%% Algorithm AS 144
if strcmp(algorithm,'all') || strcmp(algorithm,'144')
    %  Initialize vector to be permuted.
    
    %% OLD code
    %  nvect=1:ntotal;
    %  Permute vector.
    %    nnvect=nvect;
    %     ntemp = ntotal;
    %     for i = 1 : ntotal
    %         noct = floor ( rand(1,1) * ntemp + 1.0 );
    %         nvect(i) = nnvect(noct);
    %         nnvect(noct) = nnvect(ntemp);
    %         ntemp = ntemp - 1;
    %     end
    
    %% NEW CODE
    % generate a vector of pseudorandom scalar integers between 1 and ntotal
    nvect=randi(ntotal,1,ntotal);
    
    % Initialize matrix which will contain required simulated contingency
    % table using algorithm AS144
    matrix144=zeros(I,J);
    ii = 1;
    for i = 1 : I
        limit = nrowt(i);
        for k = 1 : limit
            for j = 1 : J
                if nvect(ii) <= nsubt(j)
                    ii = ii + 1;
                    matrix144(i,j) = matrix144(i,j) + 1;
                    break
                end
            end
        end
    end
    out.m144=matrix144;
end

%% Algorithm AS 159
if strcmp(algorithm,'all') || strcmp(algorithm,'159')
    % Initialize matrix which will contain required simulated contingency
    % table using algorithm AS159
    matrix159=zeros(I, J);
    
    
    % use logarithm of the gamma function to compute log factorials
    fact=gammaln(1:ntotal+1);
    
    %      % OLD inefficient implementation to calculate log-factorials.
    %    fact = zeros(ntotal,1);
    %     x = 0.0;
    %     fact(1) = 0.0;
    %     for i = 1 : ntotal
    %         x = x + log ( i );
    %         fact(i+1) = x;
    %     end
    
    %
    %  Construct a random matrix.
    %
    jwork(1:J-1) = ncolt(1:J-1);
    
    jc = ntotal;
    
    % Generate a set of (I-1)*(J-1) numbers  from uniform
    rchk=rand((I-1)*(J-1),1);
    
    for l = 1 : I - 1
        
        nrowtl = nrowt(l);
        ia = nrowtl;
        ic = jc;
        jc = jc - nrowtl;
        
        for m = 1 : J - 1
            
            id = jwork(m);
            ie = ic;
            ic = ic - id;
            ib = ie - ia;
            ii = ib - id;
            %
            %  Test for zero entries in matrix.
            %
            
            if  ie == 0
                ia = 0;
                matrix159(l,m:J) = 0;
                break
            end
            
            %  Generate a pseudo-random number from uniform distribution.
            % r=rand(1,1);
            r=rchk(l*m);
            
            %  Compute the conditional expected value of MATRIX(L,M).
            done1 = 0;
            
            while  1
                nlm = floor ( ia * id / ie + 0.5 );
                iap = ia + 1;
                idp = id + 1;
                igp = idp - nlm;
                ihp = iap - nlm;
                nlmp = nlm + 1;
                iip = ii + nlmp;
                x = exp ( fact(iap) + fact(ib+1) + fact(ic+1) + fact(idp) - ...
                    fact(ie+1) - fact(nlmp) - fact(igp) - fact(ihp) - fact(iip) );
                
                if ( r <= x )
                    break
                end
                
                sumprb = x;
                y = x;
                nll = nlm;
                lsp = 0;
                lsm = 0;
                
                %  Increment entry in row L, column M.
                while ~lsp
                    
                    j = ( id - nlm ) * ( ia - nlm );
                    
                    if  j == 0
                        
                        lsp = 1;
                        
                    else
                        
                        nlm = nlm + 1;
                        x = x * j / ( nlm * ( ii + nlm ) );
                        sumprb = sumprb + x;
                        
                        if  r <= sumprb
                            done1 = 1;
                            break
                        end
                        
                    end
                    
                    done2 = 0;
                    
                    while ( ~lsm )
                        %
                        %  Decrement the entry in row L, column M.
                        %
                        j = nll * ( ii + nll );
                        
                        if  j == 0
                            lsm = 1;
                            break
                        end
                        
                        nll = nll - 1;
                        y = y * j / ( ( id - nll ) * ( ia - nll ) );
                        sumprb = sumprb + y;
                        
                        if  r <= sumprb
                            nlm = nll;
                            done2 = 1;
                            break
                        end
                        
                        if  ~lsp
                            break
                        end
                        
                    end
                    
                    if  done2
                        break
                    end
                    
                end
                
                if  done1
                    break
                end
                
                if  done2
                    break
                end
                r=rand(1,1);
                r = sumprb * r;
                
            end
            
            matrix159(l,m) = nlm;
            ia = ia - nlm;
            jwork(m) = jwork(m) - nlm;
            
        end
        
        matrix159(l,J) = ia;
        
    end
    %  Compute the last row.
    matrix159(I,1:J-1) = jwork(1:J-1);
    matrix159(I,J) = ib - matrix159(I,J-1);
    out.m159=matrix159;
end

end

%FScategory:MULT-Categorical
