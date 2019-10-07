%% Example 1
% FSreg called using an input table and all default arguments
% In this case by default the last column of the table is used as response
% and function FSR is called
load hald
ingredients(1:5,4)=1;
ingredients(6:end-3,4)=2;
ingredients(end-2:end,4)=3;
% lm = fitlm(ingredients,heat)

% Case 1: input is just a table 
tbl=array2table([ingredients heat]);
out=FSreg(tbl);

%% Example 2
% FSreg is called using an input table and as second argument modelspec
load carsmall
MPG(isnan(MPG))=20;
d = table(Weight,Model_Year,MPG);
out=FSreg(d,'quadratic');

%% Example 3
% FSreg is called using an input table and a series of name pairs but
% omitting modelspec
% In the optional name/pairs below we define the column associated with the
% response variable
out=FSreg(tbl,'ResponseVar',3);

%%  Example 4: one of the regressors is defined categorical. 
% If one of the regressors is recognized as categorical then it
% is automatically expanded
load carsmall
MPG(isnan(MPG))=20;
Model_Year=categorical(Model_Year);
d = table(Weight,Model_Year, MPG);
out=FSreg(d,'linear');


%% Example 5: use of option CategoricalVars 
% FSreg called using an input table. Categorical variables are directly
% specified in input option CategoricalVars
load carsmall
MPG(isnan(MPG))=20;
tbl = table(Weight,Model_Year, MPG);
out=FSreg(tbl,'CategoricalVars',logical([0,1,0]));

%% Example 6: input is a table and model is supplied using a formula
load carsmall
MPG(isnan(MPG))=20;
d = table(Weight,Model_Year,MPG);
out=FSreg(d,'MPG ~ Weight + Weight^2');


%% Example 7: input are matrices and not tables
n=100;
y=randn(n,1);
X=randn(n,3);
out=FSreg(X,y,'quadratic');

%% Example 8: input are matrices and VarNames is supplied
n=100;
y=randn(n,1);
X=randn(n,3);
VarNames={'X1' 'X2' 'X3' 'myy'};
out=FSreg(y,X,'quadratic','VarNames',VarNames);


%% Example 9: input are matrices and model is supplied using a formula
n=100;
y=randn(n,1);
X=randn(n,3);
VarNames={'X1' 'Weight' 'X3' 'MPG'};
out=FSreg(X,y,'MPG ~ Weight + Weight^2','VarNames',VarNames);


%% Example 10: input are matrices model is supplied using a formula and option CategoricalVars is used 
n=100;
y=randn(n,1);
X=[randn(n,2) randi([1 3],n,1)];

VarNames={'X1'; 'Weight'; 'X3cat'; 'MPG'};
out=FSreg(X,y,'MPG ~ X3cat + Weight^2','VarNames',VarNames,'CategoricalVars',logical([0 0 1 0]));

%% Example 11: input are matrices and no additional input is supplied
n=100;
y=randn(n,1);
X=randn(n,2);
out=FSreg(X,y);

%% Example 12: input are matrices and input structure Control is used
n=100;
y=randn(n,1);
X=randn(n,2);
Control=struct;
Control.nsamp=100;
Control.init=round(n/2);
out=FSreg(X,y,'Control',Control);

%% Example of use of estimator Sregeda
Monitoring=true;
Estimator='S';
out=FSreg(X,y,'Estimator',Estimator,'Monitoring',Monitoring);


