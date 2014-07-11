%demoMixSim prepares connection with R

% Description: MixSim allows simulating mixtures of Gaussian
% distributions with different levels of overlap between mixture
% components.  Pairwise overlap, defined as a sum of two
% misclassification probabilities, measures the degree of
% interaction between components and can be readily employed to
% control the clustering complexity of datasets simulated from
% mixtures. These datasets can then be used for systematic
% performance investigation of clustering and finite mixture
% modeling algorithms. Among other capabilities of MixSim, there
% are computing the exact overlap for Gaussian mixtures,
% simulating Gaussian and non-Gaussian data, simulating outliers
% and noise variables, calculating various measures of agreement
% between two partitionings, and constructing parallel
% distribution plots for the graphical display of finite mixture
% models.
%% Test main MixSim functions

R_seed = 1234;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

fail = 1;
while fail
    Q = MixSim(5, 7, 'BarOmega' , 0.01, 'MaxOmega' , 0.05, 'R_seed', R_seed);
    fail = Q.fail;
end

[OmegaMap, BarOmega, MaxOmega, rcMax] = overlap(5, 7, Q.Pi, Q.Mu, Q.S);
disp(OmegaMap);
disp(BarOmega);
disp(MaxOmega);
disp(rcMax);

%% Example 1 of Section 3.1

R_seed = 1234;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

ex1 = MixSim(4, 5, 'BarOmega' , 0.05, 'MaxOmega' , 0.15, 'R_seed', R_seed);

%% Example 2 of Section 3.1

R_seed = 1234;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

ex2 = MixSim(3, 2, 'MaxOmega' , 0.1, 'sph' , true, 'PiLow' , 0.1, 'R_seed', R_seed);

%% Example 3 of Section 3.1

R_seed = 1234;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

ex3 = MixSim(2, 4, 'BarOmega' , 0.05, 'sph' , true, 'hom' , true,  'tol', 1e-10, 'int', [0 10], 'R_seed', R_seed);

%% Example for Section 3.2

% iris data
Y=load('ir.txt');
p=size(Y,2);
gr=repmat(1:3,50,1);
id=gr(:);
spmplot(Y,id);
K = max(id);

%estimate mixture parameters
t = tabulate(id);
Pi = t(:,3);
Mu = grpstats(Y,id,{'mean'});
S(:,:,1) = cov(Y(id==1,:));
S(:,:,2) = cov(Y(id==2,:));
S(:,:,3) = cov(Y(id==3,:));

% overlap
[OmegaMap, BarOmega, MaxOmega, rcMax] = overlap(K, p, Pi, Mu, S);

%% Example 1 of Section 3.3 plot (a)

R_seed = 1235;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q1a = MixSim(5, 2, 'MaxOmega' , 0.20, 'BarOmega' , 0.05, 'R_seed', R_seed);
[A1a , id1a] = simdataset(500, Q1a.Pi, Q1a.Mu, Q1a.S, 'R_seed', R_seed);

gscatter(A1a(:,1),A1a(:,2),id1a);

%% Example 1 of Section 3.3 plot (b)

R_seed = 1238;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q1b = MixSim(5, 2, 'MaxOmega' , 0.20, 'BarOmega' , 0.05, 'R_seed', R_seed);
[A1b , id1b] = simdataset(500, Q1b.Pi, Q1b.Mu, Q1b.S, 'R_seed', R_seed);

gscatter(A1b(:,1),A1b(:,2),id1b);

%% Example 2 of Section 3.3 plot (b)

R_seed = 1238;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q2 = MixSim(5, 2, 'MaxOmega' , 0.20, 'BarOmega' , 0.05, 'R_seed', R_seed);
[A2 , id2] = simdataset(500, Q2.Pi, Q2.Mu, Q2.S, 'R_seed', R_seed);

gscatter(A2(:,1), A2(:,2), id2);

%% Example 2 of Section 3.3 plot (c)

R_seed = 1238;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q3 = MixSim(3, 2, 'MaxOmega' , 0.1, 'int', [0.2 1], 'R_seed', R_seed);
[A3 , id3] = simdataset(300, Q3.Pi, Q3.Mu, Q3.S, 'lambda' , [0.1 10],'R_seed', R_seed);

gscatter(A3(:,1), A3(:,2), id3);

%% Example 2 of Section 3.3 plot (d)

R_seed = 1238;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q4 = MixSim(3, 2, 'MaxOmega' , 0.1, 'int', [0.2 1], 'R_seed', R_seed);
[A4 , id4] = simdataset(300, Q4.Pi, Q4.Mu, Q4.S, 'lambda' , [10 10],'R_seed', R_seed);

gscatter(A4(:,1) , A4(:,2), id4);

%% Example 3 of Section 3.3 plot (a)

R_seed = 1234;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q5 = MixSim(4, 2, 'BarOmega' , 0.01, 'R_seed', R_seed);
[A5 , id5] = simdataset(500, Q5.Pi, Q5.Mu, Q5.S, 'nout',10,'R_seed', R_seed);
gscatter(A5(:,1), A5(:,2), id5);

%% Example 3 of Section 3.3 plot (b)

R_seed = 1237;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q6 = MixSim(4, 2, 'BarOmega' , 0.01, 'R_seed', R_seed);
[A6 , id6] = simdataset(500, Q6.Pi, Q6.Mu, Q6.S, 'nout',10,'R_seed', R_seed);
gscatter(A6(:,1), A6(:,2), id6);

%% Example 4 of Section 3.3 plot (c)

R_seed = 1235;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q7 = MixSim(4, 1, 'MaxOmega' , 0.1, 'R_seed', R_seed);
% [A7 , id7] = simdataset(300, Q7.Pi, Q7.Mu, Q7.S, 'nnoise',1,'R_seed', R_seed);
% [A7 , id7] = simdataset(300, Q7.Pi, Q7.Mu, Q7.S, 'nnoise',1,'nout',20,'alpha',0.1,'int',[10 50],'lambda',[0.5 0.8],'R_seed', R_seed,'maxout',100);
% gscatter(A7(:,1), A7(:,2), id7);

%% Example 4 of Section 3.3 plot (c) bis
% Same as cell above but with all options for simdataset
R_seed = 1235;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q7 = MixSim(4, 1, 'MaxOmega' , 0.1, 'R_seed', R_seed);

[A7 , id7] = simdataset(300, Q7.Pi, Q7.Mu, Q7.S, 'nnoise',1,'nout',20,'alpha',0.1,'int',[10 50],'lambda',[0.5 0.8],'R_seed', R_seed,'maxout',100);

gscatter(A7(:,1), A7(:,2), id7);
% To run the code above in R please use the following sintax
% set.seed(1235, kind='Mersenne-Twister' , normal.kind = 'Inversion')
% Q <- MixSim(MaxOmega = 0.1, K = 4, p = 1)
% A <- simdataset(n = 300, Pi = Q$Pi, Mu = Q$Mu, S = Q$S,
%                 n.noise = 1, n.out = 20, alpha = 0.1,
%                 int = c(10,50), lambda = c(0.5,0.8), max.out = 100)


%% Example 4 of Section 3.3 plot (d)

R_seed = 1236;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q8 = MixSim(4, 1, 'MaxOmega' , 0.1, 'R_seed', R_seed);
[A8 , id8] = simdataset(300, Q8.Pi, Q8.Mu, Q8.S, 'nnoise',1,'R_seed', R_seed);
gscatter(A8(:,1), A8(:,2), id8);


%% Example of Section 3.4 plot (a)

% this cell demos the pdplot: not yet implemented

% iris data
Y=load('ir.txt');
p=size(Y,2);
gr=repmat(1:3,50,1);
id=gr(:);
spmplot(Y,id);
K = max(id);

%estimate mixture parameters
t = tabulate(id);
Pi = t(:,3);
Mu = grpstats(Y,id,{'mean'});
S(:,:,1) = cov(Y(id==1,:));
S(:,:,2) = cov(Y(id==2,:));
S(:,:,3) = cov(Y(id==3,:));

%pdplot, not yet implemented
%pdplot(Pi, Mu, S);

%% Example of Section 3.4 plot (b)
R_seed = 1234;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q9 = MixSim(6, 4, 'BarOmega' , 0.001, 'R_seed', R_seed);
%pdplot, not yet implemented
%pdplot(Q9.Pi, Q9.Mu, Q9.S);

%% Example of Section 3.4 plot (c)

R_seed = 1232;
if R_seed
    [~] = evalR(['set.seed(' num2str(R_seed) ', kind=''Mersenne-Twister'', normal.kind = ''Inversion'')']);
end

Q10 = MixSim(6, 4, 'BarOmega' , 0.05, 'R_seed', R_seed);
%pdplot, not yet implemented
%pdplot(Q10.Pi, Q10.Mu, Q10.S);


%% close connection with R

if R_seed
    closeR;
end