Y=load('hbk.txt');
[n,p]=size(Y);
% alpha = asymptotic rejection point
alpha=0.1;

medy=median(Y);
Yst=zscoreFS(Y);
% mu = initial centroid which comes from MVE
mu=[-0.08381246123673669 -0.12892364167151721 -0.13115098100933800];
% v= = intial shape matrix which comes from MVE (det(v)=1)
 v=[1.2257311771795074  0.4373758448293255  0.5945330282040835;
 0.4373758448293255  1.3145504154106444 -0.3478278166515631;
 0.5945330282040835 -0.3478278166515631  1.2632798585812688];
maxiter=120;
tol=1e-05;

outf=iterrocke(Yst, mu, v, maxiter, tol, alpha);
covv=outf.cov;
mah=outf.mah;
iter=outf.iter;
center=outf.center;
    med0 = median(mah);
    qchis5 = chi2inv(0.5,p);
    covv = (med0/qchis5)*covv;
    mah = (qchis5/med0)*mah;
devin=1.4826*mad(Y,1);
    DEV1 = diag(devin);
    center  = center*DEV1+medy;
    
    covv = DEV1*covv*DEV1;
    
    crit =det(covv);
   
    outf=struct;
    outf.cov=covv;
    outf.crit=crit;
   outf.iter=iter;
   outf.mah=mah;
   
