function [medout , varargout] = medcouple(z, mcm, wmm)
%MEDCOUPLE computes the medcouple, a robust skewness estimator.
%
%<a href="matlab: docsearchFS('medcouple')">Link to the help function</a>
%
%  Required input arguments:
%
%    z : Uni-dimensional array. The elements of the array z represent a  
%        sample from a univariate distribution. They can be missing (nan)
%        of infinite (inf) values, which are excluded from the computation
%        if present.
%
%  Optional input arguments:
%
%    mcm : Scalar. Medcouple Method. By default (mcm = 0 or not specified)
%          the medcouple is computed using a fast $O(n log n)$ algorithm.
%          mcm = 1 applies a naive $O^2(n)$ method. With mcm = 2 the
%          medcouple is replaced by a quantile-based approximation.
%          Similarily, mcm = 3 uses a octile-based approximation.
%          All methods are discussed in the references reported below. The
%          fast options are the most precise; the simplified (naive)
%          algotithm is also very fast but may be less accurate, as the
%          kernel function of Hubert and Struyf is implemented in a way
%          that does not take into account possible zeros at the
%          denominator (these cases are simply not considered in the
%          calculation). The quantile and octile approximations by
%          definition are robust skewness estimators alternatives to the
%          medcouple.
%          Example - 'mcm',2
%          Data Types - single | double
%
%    wmm : Scalar. Weighted Median Method. By default (wmm=0 or not specified)
%          the weighted median is computed using function quickselectFSw,
%          which csn compute any weighted quantile. If wmm=1 the weighted
%          median is computed using the algorithm by Sven Haase
%          https://it.mathworks.com/matlabcentral/fileexchange/23077-weighted-median
%          Option wmm=2 applies the mex file quickselectFSwmex.
%          Example - 'wmm',1
%          Data Types - single | double
%
%  Output:
%
%  medout  : The medcouple of z. Scalar. This is a robust measure of  
%            skewness for the elements in z.
%
%  Optional Output:
%
%    T     : The execution time spent for computing the weighted median.
%            Scalar. The computation of the weighted median is a critical
%            step in the medcouple and it is also computational demanding.
%            For assessment purposes it is reported optionally.
%
% More About:
%
%﻿The fast algorithm for the medcouple (Brys G., Hubert M. and Struyf A., 
% 2004) is a robust measure of skewness used to adjust the whiskers of a
% boxplot (Hubert M. and Vandervierenb E., 2008) and avoid wrong
% declaration of outliers in asymmetric univariate data. The medcouple
% replaces the classical skewness (that is based on the third moment of a
% data distribution) with a scaled median difference of the left and right
% half of the distribution. More precisely, it is defined as the median of
% the kernel function
% $h(z_i,z_j)=\frac{(z_j-median(z))-(median(z)-z_i)}{z_j-z_i}$ over all
% pairs $(z_i,z_j)$ where $z_i<=median(z)$ and $z_j >= median(z)$. The MATLAB
% toolbox LIBRA (http://wis.kuleuven.be/stat/robust.html) and the R package
% robustbase (https://cran.r-project.org/web/packages/robustbase) implement
% the fast medcouple with wrapper functions to the same compiled C source,
% to maximize speed. This medcouple.m function replicates faithfully the C
% function mlmc.c in pure MATLAB code, but replaces the key computations of
% the median and weighted median (which are intensively used) with calls to
% the FSDA functions quickselectFS and quickselectFSw or, optionally, to
% the $n\log n$ solution of Sven Haase (2022) weightedMedian.m. The result
% is a readable MATLAB function that does not ﻿sacrifices performance. The
% FSDA implementation also optimizes the computation of the medcouple
% kernel function, to take into account the case where there are multiple
% values equal to the median, the denominator becomes 0 and the kernel is
% therefore defined differently: we computed it efficiently in the
% sub-function calwork, using the signum function. The consistency of the
% MATLAB and original C implementations have been carefully verified with
% simulations that ensure the same random numbers generation in the two
% environments. The code that has been used can be found in the FSDA github
% project.
%
% See also: quickselectFS, quickselectFSw
%
%
% References:
%
% Brys G., Hubert M. and Struyf A. (2004), A Robust Measure of Skewness,
% "Journal of Computational and Graphical Statistics", Vol. 13(4), pp.
% 996-1017.
%
% Hubert M. and Vandervierenb E. (2008), An adjusted boxplot for skewed
% distributions, "Computational Statistics and Data Analysis", Vol. 52, pp.
% 5186-5201.
%
% Johnson D.B. and Mizoguchi T. (1978), Selecting the Kth Element in X + Y
% and X1 + X2 + ... + Xm, "SIAM Journal of Computing", Vol. 7, pp. 147-153.
%
% Hinkley D.V. (1975), On Power Transformations to Symmetry, "Biometrika",
% Vol. 62, pp.101-111.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('medcouple')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2019-11-20 18:33:09 #$: Date of the last commit
%
%
% Examples:
%
%{
    %% the fast medcouple (default) on chi2 data.
    rng(5);
    data = chi2rnd(5,1000,1);
    mc1 = medcouple(data);
    stddev = std(data);
    figure;
    histogram(data);
    title(['Medcouple = ' num2str(mc1) ' -- Standard deviation = ' num2str(stddev)],'FontSize',16);
%}

%{
    % the naive medcouple on chi2 data
    rng(5);
    method = 1;
    mc2 = medcouple(chi2rnd(5,1000,1),method);
%}

%{ 
    %% The medcouple on small datasets.
    % For small n, it is better to compute the the medcouple also on 
    % the reflected sample -z and return (medcouple(z) - medcouple(-z))/2;
    % To work on positive values we do this with:

    z = chi2rnd(5,25,1);
    mcm = 1; % use naive method (n is small ...)
    medoutR = 0.5*(-medcouple(repmat(max(z),size(z,1),1)-z , mcm) + medcouple(z,mcm));
    medout  = medcouple(z,mcm);
    disp('Application to a non problematic dataset:');
    disp(['result with reflected option = ' num2str(medoutR) ' --- result on z only = ' num2str(medoutR)]);

%} 

%{
    %% Boxplot using medcouple to adjust for skewness.
    % Apply to symmetric and skewed data

    rng(5);

    %orientation = 'horizontal'; 
    orientation = 'vertical'; 

    % normal data with an outlier
    datain1 = 15+randn(200,1);     
    datain1 = [datain1 ; 35];      
 
    % log-normal data with an outlier
    datain2 = lognrnd(0,1,200,1);  
    datain2 = [datain2 ; 35]; 

    k  = 1.5;

    % limits on normal data + one outlier
    MC_1  = medcouple(datain1);
    Q1_1  = quantile(datain1,0.25);
    Q3_1  = quantile(datain1,0.75);
    IQR_1 = Q3_1 - Q1_1;
    DataLimLower_1 = Q1_1 - k * exp(-3.5 * MC_1) * IQR_1;
    DataLimUpper_1 = Q3_1 + k * exp(4 * MC_1)    * IQR_1;

    % limits on log-normal data + one outlier
    MC_2  = medcouple(datain2);
    Q1_2  = quantile(datain2,0.25);
    Q3_2  = quantile(datain2,0.75);
    IQR_2 = Q3_2 - Q1_2;
    DataLimLower_2 = Q1_2 - k * exp(-3.5 * MC_2) * IQR_2;
    DataLimUpper_2 = Q3_2 + k * exp(4 * MC_2)    * IQR_2;

    % boxplot
    boxplot([datain1,datain2],'Whisker',k,'Orientation',orientation,...
        'OutlierSize',10,'Symbol','*r', ...
        'labels',{'Normal(15,1)+outlier','Log-normal(0,1)+outlier'});

    title({'Boxplot adjusted for skewness' , 'with medcouple (red limits)'}, 'FontSize' , 20)
    set(gca,'FontSize',16);
    
    switch orientation
        case 'vertical'
            % add the limits (for vertical boxplot)
            xx  = xlim(gca);
            xx2 = xx/1.1;   xx2(1) = xx(2)/1.5;
            xx1 = xx/2; xx1(1) = xx(2)/3.5;
            line(xx1,[DataLimLower_1 , DataLimLower_1] , 'Color','red');
            line(xx1,[DataLimUpper_1 , DataLimUpper_1] , 'Color','red');
            line(xx2,[DataLimLower_2 , DataLimLower_2] , 'Color','red');
            line(xx2,[DataLimUpper_2 , DataLimUpper_2] , 'Color','red');
        case 'horizontal'
            % add the limits (for horizontal boxplot)
            yy  = ylim(gca);
            yy2 = yy/1.1;   yy2(1) = yy(2)/1.5;
            yy1 = yy/2; yy1(1) = yy(2)/3.5;
            line([DataLimLower_1 , DataLimLower_1] , yy1,'Color','red');
            line([DataLimUpper_1 , DataLimUpper_1] , yy1,'Color','red');
            line([DataLimLower_2 , DataLimLower_2] , yy2,'Color','red');
            line([DataLimUpper_2 , DataLimUpper_2] , yy2,'Color','red');
     end
%}

%{
    %% Assess the perfomance of the various medcouple solutions.

    % Uses the time execution of the weighted median returned in varargout

    % Tm=matlab-fast-quickselectFSw; Tw=matlab-fast-whimed
    % Ts=matlab-simplified; To=matlab-octile; Tq=matlab-quantiles;
    tm=0; tmmex=0; tw=0; ts=0; to = 0; tq=0;

    N      = 800;   % sample size
    cycles = 1000;  % number of repetitions

    for c = 1:cycles

        if mod(c,floor(cycles/10))==0, disp(['repetition n. ' num2str(c)]); end

        % Generate different random sequences (seed 896 is a good one to
        % compare consistency with the original C code on a difficult sequence)
        myseed = randi(1000,1,1);
        rng(myseed , 'twister');
        datain = rand(N,1);
        % REMARK: to generate N random integers (for example between 1 and
        % 1000) using the Mersenne Twister in both C and MATLAB, note that the
        % line below is equivalent to randi(1000,N,1):
        %datain = ceil(datain*1000);

        % Line below is to test with skewed data from the lognormal
        %datain = lognrnd(0,1,N,1);

        % just to ensure that data are stored in a column vector
        datain = datain(:);

        % $n\log n$, using quickselectFSw
        tm0 = tic;
        [MCm , t_diva] = medcouple(datain,0,0);
        tm  = tm+toc(tm0);

        % $n\log n$, using quickselectFSwmex
        tm0mex = tic;
        [MCmmex , t_diva_mex] = medcouple(datain,0,2);
        tmmex  = tmmex+toc(tm0mex);

        % $n\log n$, using Haase's weighted median (which is based on MATLAB sortrows)
        tw0 = tic;
        [MCw , t_Haase] = medcouple(datain,0,1);
        tw  = tw+toc(tw0);

        % $n^2$ "naive" algorithm, simplified
        ts0 = tic;
        MCs = medcouple(datain,1);
        ts  = ts+toc(ts0);

        % Quantiles-based approximation
        to0 = tic;
        MCo = medcouple(datain,2);
        to  = to+toc(to0);

        % Octiles-based approximation
        tq0 = tic;
        MCq = medcouple(datain,3);
        tq  = tq+toc(tq0);

    end

    disp(' ');
    disp(['OVERALL TIME EXECUTION IN ' num2str(cycles) ' MC COMPUTATION']);
    disp(['time using WM by Haase       = ' num2str(tw)]);
    disp(['time using quickselectFSw    = ' num2str(tm)]);
    disp(['time using quickselectFSwmex = ' num2str(tmmex)]);
    disp(['time using the naive         = ' num2str(ts)]);
    disp(['time using quantiles         = ' num2str(tq)]);
    disp(['time using octiles           = ' num2str(to)]);
    disp(' ');
    disp('TIME SPENT IN COMPUTING WEIGHTED AVERAGES IN A SINGLE MC COMPUTATION')
    disp(['t_Haase             = ' num2str(t_Haase)]);
    disp(['t_quickselectFSw    = ' num2str(t_diva)]);
    disp(['t_quickselectFSwmex = ' num2str(t_diva_mex)]);
%}

%% preliminaries

% Check for proper number of arguments */
if nargin > 3
    error('Wrong number of input arguments');
end

if nargin < 3 || ~ismember(wmm , [0,1,2]) || isempty(wmm)
    wmm = 0;  % default is the whimed function
end

if nargin < 2 || ~ismember(mcm , [0,1,2,3]) || isempty(mcm)
    mcm = 0;  % default is the fast method
end

% get rid of inf and nan
z(sum(~isfinite(z),2)>0,:)=[];

[n,m] = size(z);

if (m~=1 && n==1)
    error('Input must be a columnvector');
end

if (n>50000)
    warning('Array too large: you can run into memory/time problems.');
end

switch mcm
    
    case 0
        %% fast algorithm O(n log n)
        % Uses result of Johnson D.B. and Mizoguchi T. (1978)
        % Uses quickselectFS for median and quickselectFSw for weighted median
        
        % used to estimate time execution of the weighted median
        t_Haase=0; t_diva=0;
        
        % initializations
        y = zeros(n,1);
        y1 = y; y2 = y;
        %left = y; right = y;q= y; p = y; work = y; weight = y;
        epsilon = 10^(-13);
        
        % More convenient to go in decreasing order: then, use opposite numbers.
        x=-z;
        
        % compute initial ingredients, as for the naive medcouple
        
        % median
        xmed = quickselectFS(x,floor(n/2)+1);
        if mod(n,2) == 0
            xmed2 = quickselectFS(x,floor(n/2));
            xmed  = (xmed+xmed2)/2;
        end
        
        % Centred and rescaled vectors
        x = x - xmed;
        
        % Sorting: O(n log n) time
        y = sort(x);
        
        % Begin Kth pair algorithm (Johnson & Mizoguchi)
        
        if -y(1) > y(n)
            yden = -2*y(1);
        else
            yden =  2*y(n);
        end
        y = -y/yden;
        
        j=1;
        while y(j) > epsilon
            y1(j)=y(j);
            j=j+1;
        end
        
        i=1;
        while y(j) > -epsilon
            y1(j)=y(j);
            y2(i)=y(j);
            j=j+1;
            i=i+1;
        end
        
        %{
        % The two loops above are not slower than this loopless block

        jposeps = y > epsilon;
        y1(jposeps) = y(jposeps);

        jnegeps = y > -epsilon;
        jnegeps = and(~jposeps,jnegeps);

        y1(jnegeps) = y(jnegeps);
        y2(1:sum(jnegeps)) = y(jnegeps);

        j = sum(jposeps)+sum(jnegeps)+1;
        i = sum(jnegeps)+1;
        %}
        
        h1 = j-1;
        while j < n+1
            y2(i)=y(j);
            j=j+1;
            i=i+1;
        end
        
        % The initial left and right boundaries
        h2 = i-1;
        
        left   = zeros(1,h2);
        right  = left;
        q      = left;
        p      = left;
        work   = left;
        weight = left;
        
        for i=1:h2
            left(i)  = 1;
            right(i) = h1;
        end
        
        nl   = 0;     % number of entries to the left of the left boundary
        nr   = h1*h2; % number of entries to the left of the right boundary
        knew = floor(nr/2)+1; % the medcouple index
        
        IsFound=false;
        while (nr-nl>n) && ~IsFound
            % Iterate while the number of entries between the boundaries is
            % greater than the number of rows in the matrix
            
            % Compute row medians (work) and their associated weights (weight)
            j=1;
            for i=1:h2
                if left(i)<=right(i)
                    weight(j)=right(i)-left(i)+1;
                    k = left(i)+floor(weight(j)/2);
                    % work(j) = calwork(y1(k),y2(i),k,i,h1+1);
                    % lines below are faster than call to calwork
                    work(j) = (y1(k)+y2(i))/(y1(k)-y2(i));  % Kernel function
                    if isnan(work(j))                       % Kernel function
                        work(j) = -sign(k+i-h1-1);          % Kernel function
                    end
                    
                    j=j+1;
                end
            end
            
            % Approximate medcouple with the weighted median of the row medians
            
            % Weights normalised to one
            weight_jm1   = weight(1:j-1) / sum(weight(1:j-1));
            % Data
            work_jm1     = work(1:j-1);
            if wmm == 0
                % Use the weighted median implementation of FSDA. 
                % It is a O(n) algorithm.
                tt        = tic();
                trial     = quickselectFSw(work_jm1,weight_jm1,0.5);
                t_diva    = t_diva+toc(tt);
            elseif wmm == 1
                % Use the weighted median implementation by Sven Haase
                % It is a O(n log(n)) algorithm (uses sortrows).
                t           = tic();
                [~ , trial] = weightedMedianHaase(work_jm1,weight_jm1);
                t_Haase    = t_Haase+toc(t);
            else
                % quickselectFSw mex ... Optimized
                ttt = tic;
                work_copy = work_jm1;   work_copy(end+1)=999; work_copy(end)=[]; %#ok<AGROW>
                weig_copy = weight_jm1; weig_copy(end+1)=999; weig_copy(end)=[]; %#ok<AGROW>
                trial = quickselectFSwmex(work_copy,weig_copy,0.5,j-1);
                t_diva    = t_diva+toc(ttt);
            end
            
            % New tentative right boundary
            j=1;
            for i=h2:-1:1
                % while ((j<=h1) && (calwork(y1(j),y2(i),j,i,h1+1)>trial))
                %    j=j+1;
                % end
                % block below is much faster than call to calwork in while
                check = true;
                while (j<=h1 && check)
                    cwork = (y1(j)+y2(i))/(y1(j)-y2(i));    % Kernel function
                    if isnan(cwork)                         % Kernel function
                        cwork = -sign(j+i-h1-1);            % Kernel function
                    end
                    check=cwork>trial;
                    if check
                        j=j+1;
                    end
                end
                % end of block
                p(i)=j-1;
            end
            
            % New tentative left boundary
            j=h1;
            for i=1:h2
                % while ((j>=1) && (calwork(y1(j),y2(i),j,i,h1+1)<trial))
                %    j=j-1;
                % end
                % block below is much faster than call to calwork in while
                check = true;
                while (j>=1 && check)
                    cwork = (y1(j)+y2(i))/(y1(j)-y2(i));    % Kernel function
                    if isnan(cwork)                         % Kernel function
                        cwork = -sign(j+i-h1-1);            % Kernel function
                    end
                    check=cwork<trial;
                    if check
                        j=j-1;
                    end
                end
                % end of block
                q(i)=j+1;
            end
            
            % Determine which entries to discard, or if we've found the medcouple
            sump=0;
            sumq=0;
            for i=1:h2
                sump=sump+p(i);
                sumq=sumq+q(i)-1;
            end
            
            if (knew<=sump)
                for i=1:h2
                    right(i)=p(i);
                end
                nr=sump;
            else
                if (knew>sumq)
                    for i=1:h2
                        left(i)=q(i);
                    end
                    nl=sumq;
                else
                    medout=trial;
                    IsFound=true;
                    % Found the medcouple
                end
            end
        end %end while-lus*/
        
        % Did not find the medcouple, but there are very few tentative entries
        % remaining
        if ~IsFound
            j=1;
            for i=1:h2
                if (left(i)<=right(i))
                    for jj=left(i):right(i)
                        % work(j) = -calwork(y1(jj),y2(i),jj,i,h1+1);
                        % line below are faster than call to calwork
                        work(j) = -(y1(jj)+y2(i))/(y1(jj)-y2(i));   % Kernel function
                        if isnan(work(j))                           % Kernel function
                            work(j) = +sign(jj+i-h1-1);             % Kernel function
                        end
                        j=j+1;
                    end
                end
            end
            % Select the medcouple by rank (knew-nl) amongst the remaining
            % entries (the first j-1 elements of work only)
            medout = quickselectFS(work(1:j-1),knew-nl);
            medout = -medout;
        end
        
        % Store in varargout the time spent on the weighted median
        % computation, which is the expensive part of the algorithm. It is
        % provided for assessment purposes.
        if nargout==2
            if wmm==1
                varargout={t_Haase};
            else
                varargout={t_diva};
            end
        end
        
    case 1
        %% naive algorithm O(n^2)
        % Evaluates the kernel function $h$ for each couple of points
        % (xi,xj) with xi ≤ median and xj ≥ median. Therefore this
        % algorithm needs $O(n^2)$ time.
        
        x = z(~isnan(z)); % remove nan
        
        m = median(x);    % get median
        
        y = sort(x);      % sort data
        
        Ygi = y(y>=m);    % All values >= to median
        Ylj = y(y<=m);    % All values <= to median
        
        % build all possible combinations
        Yi = repmat(Ygi ,1, length(Ylj)); 
        Yj = repmat(Ylj',length(Ygi), 1); 
        
        % formula 2.2 of Brys G., Hubert M. and Struyf A. (2004)
        h  = (((Yi-m)-(m-Yj))./(Yi-Yj));
        
        % the case of multiple values tied to the median is not addressed
        medout = median(h(:),'omitnan');
        
    case 2
        %% approximation based on octile skewness
        %  Hinkley D.V. (1975)
        
        Q = quantile(z,[0.125 0.50 0.875]);
        Q125 = Q(1); Q50 = Q(2); Q875 = Q(3);
        medout = ((Q875 - Q50) - (Q50 - Q125)) / (Q875 - Q125);
        
    case 3
        %% approximation based on quartile skewness
        %  Hinkley D.V. (1975)
        
        Q = quantile(z,[0.25 0.50 0.75]);
        Q25 = Q(1); Q50 = Q(2); Q75 = Q(3);
        medout = ((Q75 - Q50) - (Q50 - Q25)) / (Q75 - Q25);
        
end

end

% ---------------------------------------------------------
%% internal functions
% ---------------------------------------------------------

function [wdSort , wMed , j] = weightedMedianHaase(D,W)

% Function to compute the weighted median of n elements in $O(n \log n)$
% worst-case time. It uses the sorting method implemented by MATLAB
% function sortrows. The code, by Sven Haase, is distributed in the
% Mathworks community website.
%
% Reference:
% Sven Haase (2022). Weighted median. MATLAB Central File Exchange.
% Retrieved January 30, 2022.
% (https://www.mathworks.com/matlabcentral/fileexchange/23077-weighted-median)
%
% Algorithm:
%
% For n numbers x_1,...,x_n with positive weights w_1,...,w_n,
% (sum of all weights equal to one) the weighted median is defined as
% the element x_k, such that:
%           --                        --
%           )   w_i  <= 1/2   and     )   w_i <= 1/2
%           --                        --
%        x_i < x_k                 x_i > x_k
%
%
% Input:    D ... matrix of observed values
%           W ... matrix of weights, W = ( w_ij )
% Output:   wMed ... weighted median
%
%
% Copyright (c) 2009, Sven Haase
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if nargin ~= 2
    error('weightedMedian:wrongNumberOfArguments', ...
        'Wrong number of arguments.');
end

if size(D) ~= size(W)
    error('weightedMedian:wrongMatrixDimension', ...
        'The dimensions of the input-matrices must match.');
end

% normalize the weights, such that: sum ( w_ij ) = 1
% (sum of all weights equal to one)

WSum = sum(W(:));
W = W / WSum;

% (line by line) transformation of the input-matrices to line-vectors
d = reshape(D',1,[]);
w = reshape(W',1,[]);

% sort the vectors
A = [d' w'];
ASort = sortrows(A,1);

dSort = ASort(:,1)';
wSort = ASort(:,2)';

% vector for cumulative sums of the weights
sumVec = zeros(length(wSort),1); 
for i = 1:length(wSort)
    sumVec(i) = sum(wSort(1:i));
end

wMed = [];
j = 0;

while isempty(wMed)
    j = j + 1;
    if sumVec(j) >= 0.5
        wMed = dSort(j);    % value of the weighted median
    end
end

wdSort = [dSort ; wSort];

% final test to exclude errors in calculation
if ( sum(wSort(1:j-1)) > 0.5 ) && ( sum(wSort(j+1:length(wSort))) > 0.5 )
    error('weightedMedian:unknownError', ...
        'The weighted median could not be calculated.');
end
end

%% ----------------------------------------------------------------------
%  An unused function, reported here as a reference as its role in the
%  medcouple computation is important.
%  ----------------------------------------------------------------------
%
% calwork(a,b,ai,bi,ab) returns the kernel values needed to compute the mc.
% This version optimises the original one ported from C, which follows
% below. It is unused because this calculation is repeated many times and
% for optimizating the execution it is preferable to report these three
% lines in the body of the medcouple code.
%
% The optimization is obtained as follows. The kernel function is
% $h(z_i,z_j)=\frac{(z_j-median(z))-(median(z)-z_i)}{z_j-z_i}$ and the
% medcouple is the median of the $h(z_i,z_j)$ values computed over all
% pairs $(z_i,z_j)$ where $z_i<=median(z)$ and $z_j >= median(z)$. Function
% calwork_oric below implements the original implementation by the
% medcouple authors in C. To address the singularity occurring at the
% denominator where there are multiple values tied to the median (equation
% 2.3 of Brys et al 2004), we have re-defined the function by exploiting
% the signum function of MATLAB.

function cwork = calwork(a,b,ai,bi,ab) %#ok<DEFNU>
cwork = (a+b)/(a-b);
if isnan(cwork)
    cwork = -sign(ai+bi-ab);
end

end

%% ----------------------------------------------------------------------

% calwork_oric(a,b,ai,bi,ab) is the original version of calwork, ported
% from C without major changes
function cwork = calwork_oric(a,b,ai,bi,ab) %#ok<DEFNU>

if abs(a-b) < 2.0*10^(-13)
    if (ai+bi == ab)
        cwork = 0;
    else
        if ai+bi < ab
            cwork = 1;
        else
            cwork = -1;
        end
    end
else
    cwork = (a+b)/(a-b);
end

end

%FScategory:UTISTAT


