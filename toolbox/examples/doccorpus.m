function [docs,docname,theme,themename]=doccorpus(varargin)
%doccorpus collects the help pages of FSDA and MATLAB used to illustrate the clustering of document embeddings
%
%<a href="matlab: docsearchFS('doccorpus')">Link to the help function</a>
%
% Required input arguments:
%
%  No required input arguments.
%
% Optional input arguments:
%
%   maxdoc: maximum number of documents for each MATLAB theme. Scalar.
%           When a theme has more functions, an evenly spaced sample is
%           taken, so that the themes have a similar size.
%               Example - 'maxdoc',30
%               Data Types - double
%
% Output:
%
%   docs :     documents. Cell array of n-by-1 character vectors.
%              Help text of each function.
%   docname :  names of the functions. Cell array of n-by-1 character vectors.
%   theme :    theme of each document. Vector of length n with integers
%              from 1 to numel(themename).
%   themename: names of the themes. Cell array of character vectors.
%              Themes 1-4 are categories of FSDA (regression, robust
%              clustering, distributions, utilities) and the others are
%              themes of MATLAB with no relation to robust statistics.
%
% More About:
%
% The documents are the help pages of the functions which are installed:
% FSDA functions are selected using the category written in the last line
% of each file (%FScategory:), and the text before the code and before the
% examples is used. The functions of MATLAB are taken from lists of names
% which are in this file, using the command help. Names which are not
% installed (for example because a toolbox is missing) are skipped, so that
% the number of documents depends on the toolboxes of the user. The lists
% can be modified or extended with new themes.
%
% See also: docembed, tclustvmf, examples_directional
%
% References:
%
% Salton G., Wong A., Yang C.S. (1975), A vector space model for automatic
% indexing, Communications of the ACM, 18, 613-620.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('doccorpus')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Documents and their themes.
    [docs,docname,theme,themename]=doccorpus;
    disp(array2table(accumarray(theme,1)','VariableNames',strrep(themename,' ','_')))
%}

%% Beginning of code
maxm=40;
if nargin>0
    if nargin~=2 || ~strcmp(varargin{1},'maxdoc')
        error('FSDA:doccorpus:WrongOption','The only option is: ''maxdoc'',value')
    end
    maxm=varargin{2};
end

% 1. Categories of FSDA
fsdacat={'REG-Regression','CLUS-RobClaMULT','ProbDist','UTIGEN'};
maxfsda=[35 35 35 15];
root=fileparts(fileparts(which('tclust')));       % folder of the FSDA toolbox
files=dir(fullfile(root,'**','*.m'));
skip={'+aux','examples','FSDAdemos','apps','utilities_help','doc'};
docs={}; theme=[]; docname={};
count=zeros(1,numel(fsdacat));
for i=1:numel(files)
    % part of the path below the FSDA folder, to skip some subfolders
    rel=[files(i).folder(numel(root)+1:end) filesep];
    if any(contains(rel,strcat(filesep,skip,filesep)))
        continue
    end
    txt=fileread(fullfile(files(i).folder,files(i).name));
    cc=regexp(txt,'%FScategory:(\S+)','tokens','once');
    if isempty(cc)
        continue
    end
    c=find(strcmp(cc{1},fsdacat),1);
    if isempty(c) || count(c)>=maxfsda(c)
        continue
    end
    % help text: before the code and before the examples
    ib=strfind(txt,'%% Beginning of code');
    if isempty(ib), continue, end
    txt=txt(1:ib(1)-1);
    ie=strfind(txt,'% Examples');
    if ~isempty(ie), txt=txt(1:ie(1)-1); end
    count(c)=count(c)+1;
    docs{end+1,1}=txt;                                  %#ok<AGROW>
    docname{end+1,1}=files(i).name(1:end-2);             %#ok<AGROW>
    theme(end+1,1)=c;                                    %#ok<AGROW>
end

% 2. Themes of MATLAB
matlabtheme={'strings','linear algebra','files','image processing', ...
    'signal processing','optimization'};
matlabfun={ ...
    {'strsplit','strjoin','strtrim','strrep','strfind','strcat','regexprep','regexp', ...
     'upper','lower','strcmp','strcmpi','strncmp','strncmpi','deblank','blanks', ...
     'num2str','str2double','int2str','sprintf','sscanf','strtok','contains','startsWith','endsWith'}, ...
    {'eig','svd','qr','chol','lu','inv','pinv','expm','logm','sqrtm','kron','null','orth', ...
     'cond','rank','det','trace','norm','schur','hess','balance','cross','dot','rcond','linsolve'}, ...
    {'fopen','fclose','fread','fwrite','fgetl','fgets','fprintf','fscanf','textscan', ...
     'readtable','writetable','csvread','csvwrite','dlmread','dlmwrite','importdata', ...
     'load','save','fseek','ftell','feof','frewind','readmatrix','writematrix','fileread'}, ...
    {'imread','imwrite','imshow','imshowpair','imfinfo','im2double','im2single','im2uint8', ...
     'im2uint16','im2int16','im2bw','imbinarize','graythresh','adaptthresh','rgb2gray', ...
     'im2gray','ind2rgb','rgb2ind','hsv2rgb','rgb2hsv','lab2rgb','rgb2lab','makecform', ...
     'applycform','imadjust','imadjustn','histeq','adapthisteq','imhist','imhistmatch', ...
     'imhistmatchn','stretchlim','imnoise','imsharpen','imreducehaze','imflatfield', ...
     'imlocalbrighten','localcontrast','locallapfilt','localtonemap','decorrstretch', ...
     'imfilter','imgaussfilt','imgaussfilt3','medfilt2','medfilt3','wiener2','nlfilter', ...
     'roifilt2','ordfilt2','entropyfilt','rangefilt','stdfilt','imboxfilt','imbilatfilt', ...
     'imguidedfilter','imgaborfilt','fspecial','padarray','imresize','imrotate', ...
     'imtranslate','imwarp','imref2d','imref3d','affine2d','affine3d','projective2d', ...
     'fitgeotrans','imregtform','imregister','imregdemons','imfuse','imgradient', ...
     'imgradientxy','edge','imfill','imclearborder','bwareaopen','bwareafilt','bwconncomp', ...
     'bwlabel','bwlabeln','bwperim','bwboundaries','bwmorph','bwmorph3','regionprops', ...
     'regionprops3','bwpropfilt','bwconvhull','bwdist','bwdistgeodesic','watershed', ...
     'imextendedmin','imextendedmax','imhmin','imhmax','imregionalmin','imregionalmax'}, ...
    {'findpeaks','findchangepts','finddelay','alignsignals','xcorr','xcov','cpsd', ...
     'mscohere','tfestimate','periodogram','pwelch','pmtm','plomb','pspectrum', ...
     'spectrogram','stft','istft','tfridge','wvd','fsst','ifsst','cwt','icwt','wt','fwt', ...
     'swt','modwt','imodwt','modwtmra','wfilters','dwt','idwt','wavedec','waverec', ...
     'wavedec2','waverec2','dwt2','idwt2','detcoef','appcoef','wrcoef','upcoef', ...
     'downcoef','wavedec3','waverec3','dwt3','idwt3','db2mag','db2pow','mag2db', ...
     'pow2db','freqz','grpdelay','phasez','impz','stepz','zplane','designfilt','butter', ...
     'cheby1','cheby2','ellip','besself','besselap','buttord','cheb1ord','cheb2ord', ...
     'ellipord','fir1','fir2','firls','firpm','firpmord','kaiserord','remezord', ...
     'hamming','hann','hanning','blackman','blackmanharris','bartlett','triang', ...
     'rectwin','kaiser','tukeywin','flattopwin','gausswin','nuttallwin','parzenwin', ...
     'slepian','dpss','resample','decimate','interp','upfirdn','filtfilt','filter','sosfilt'}, ...
    {'fmincon','fminunc','fminsearch','fminbnd','fzero','fsolve','lsqcurvefit', ...
     'lsqnonlin','lsqlin','lsqnonneg','linprog','intlinprog','quadprog','coneprog', ...
     'fgoalattain','fminimax','fseminf','optimoptions','optimset','optimproblem', ...
     'optimvar','optimexpr','solve','prob2struct','optimget','optimplotfval', ...
     'optimplotx','optimplotfirstorderopt','optimplotfunccount', ...
     'optimplotconstrviolation','optimplotstepsize','optimplotresnorm', ...
     'optimplotlsqfit','optimplotmilp','ga','gamultiobj','patternsearch', ...
     'particleswarm','simulannealbnd','surrogateopt','paretosearch','GlobalSearch', ...
     'MultiStart','createOptimProblem','GlobalOptimSolution','CustomStartPointSet', ...
     'RandomStartPointSet','SobolSet','LatinHypercubeSet','GlobalSearchOptions', ...
     'MultiStartOptions','optimizableVariable','bayesopt','fitrgp','lhsdesign', ...
     'lhsnorm','sobolset','haltonset'}};
for t=1:numel(matlabfun)
    txts={}; nm={};
    for i=1:numel(matlabfun{t})
        txt=help(matlabfun{t}{i});
        % names which are not installed give an empty or a very short text
        if numel(strtrim(txt))>=60
            txts{end+1}=txt;                              %#ok<AGROW>
            nm{end+1}=matlabfun{t}{i};                    %#ok<AGROW>
        end
    end
    m=numel(nm);
    if m>maxm
        s=round(linspace(1,m,maxm));
        txts=txts(s); nm=nm(s);
    end
    docs=[docs; txts(:)];                                 %#ok<AGROW>
    docname=[docname; nm(:)];                             %#ok<AGROW>
    theme=[theme; repmat(numel(fsdacat)+t,numel(nm),1)];  %#ok<AGROW>
end
themename=[{'FSDA regression','FSDA robust clustering','FSDA distributions','FSDA utilities'} ...
    matlabtheme];
end
%FScategory:CLUS-Sphere
