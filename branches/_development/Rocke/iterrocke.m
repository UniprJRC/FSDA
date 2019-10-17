function outf=iterrocke(x, center, cov, maxiter, tol, alpha)
[n,p]=size(x);
delta = 0.5*(1-p/n);

% disesca returns the estimate of the scale
disres=disesca(x,center,cov,alpha,delta);
sigmin=disres.sig;
dismin=disres.dis;
mumin=center;
vmin=cov;
it =1;
decre = tol+1;
while (decre > tol && it <= maxiter)
    
    % do reweighting (mumin and vmin are the initial values, dismin =initial vector of MD)
    des1 = descen1(x, mumin, vmin, alpha, sigmin, dismin);
    mu1 =des1.mui;
    v1 = des1.vi;
    dis1 = des1.dis;
    sig1 = des1.sig;
    decre1 = (sigmin-sig1)/sigmin;
    
    
    if (decre1 <= tol)
        
        grad =gradien(x, mu1, v1, mumin, vmin, sigmin, dismin, alpha);
        mu2 =grad.mu;
        v2 = grad.v;
        sig2 = grad.sig;
        dis2=grad.dis;
        decre2 = (sigmin-sig2)/sigmin;
        decre = max(decre1,decre2);
        if (decre1 >= decre2)
            mumin = mu1;
            vmin = v1;
            dismin = dis1;
            sigmin = sig1;
        else
            mumin = mu2;
            vmin = v2;
            dismin = dis2;
            sigmin = sig2;
        end
    else
        decre = decre1;
        mumin = mu1;
        vmin = v1;
        dismin = dis1;
        sigmin = sig1;
    end
    it = it+1;
end
outf=struct;
outf.center=mumin;
outf.cov=vmin;
outf.scale=sigmin;
outf.mah=dismin;
outf.iter=it;

    function out=disesca(x,mu,v,alpha,delta)
        % returns in structure out the estimate of the scale
        % 1) Compute MD from using initial mu and v (shape matrix)
        % rescale MD multiplying by |v|^(1/p)
        % return initial v multiplied by |v|^(-1/p)
        n =size(x,1);
        p =size(x,2);
        p2 = (1/p);
        p3 = -p2;
        % dis11 = Mahalanobis distances from shape matrix
        dis11 =mahalFS(x,mu,v);
        detv=det(v);
        
        dis11 = dis11*(detv^p2);
        v11 = v*(detv^p3);
            
        % sig is the new estimate of the scale
        % using rescaled MD
        sig = escala(p, dis11, alpha, delta);
        
        out.dis=dis11;
        out.sig=sig;
        out.v=v11;
        out.det=detv;
    end

    function sig1=escala(p,dis,alpha,delta)
        % Initial interval for scale estimate
        % x0 will contain smaller and upper values of the interval to search
        % for the solution of the equation
        % In other words the equation which has to be solved (the unknown
        % value is sig)
        % mean(rho(z,p,alpha))=delta
        % where z =dis/sig and dis is the vector of Mahalanobis distances
        % computed from a shape matrix (that is a matrix whose determinant
        % is 1)
             
        
        sig=median(dis);
        sigini=sig;
        a=sig/100;
        ff=fespro(a,p,dis,alpha,delta);
        while ff<0
            a=a/10;
            ff=fespro(a,p,dis,alpha,delta);
        end
        
        b = sig*10;
        ff=fespro(b,p,dis,alpha,delta);
        while ff>0
            b=b*10;
            ff=fespro(b,p,dis,alpha,delta);
        end
        x0=[a b];
        
        fun = @(sig) fespro(sig,p,dis,alpha,delta);
        
        sig1 = fzero(fun,x0);
        
        
        
        maxiterTMP=100;
        tolTMP=1e-12;
        sc = sigini;
loop = 0;
err = 1;
kc=delta;
while  (( loop < maxiterTMP ) && (err > tolTMP))
    % scale step: see equation 7.119 of Huber and Ronchetti, p. 176
    % scalenew = scaleold *(1/n)*\sum  \rho(u_i/scaleold) / kc
    
%             z = dis/sig;
%         sfun = mean(rhork2(z, p, alpha)) - delta;

    % scnew = sc*sqrt( mean(TBrho(u/sc,c)) / kc);
    scnew = sc*sqrt( mean(rhork2(dis/sc,p,alpha)) / kc);
    
    % Note that when there is convergence 
    % sqrt( mean(TBrho(u/sc,c)) / kc) tends to 1 (from below)
    % disp([loop sc sqrt( mean(TBrho(u/sc,c)) / kc)])
    
    err = abs(scnew/sc - 1);
    sc = scnew;
    % disp(sc)
    loop = loop+1;
end

    end




    function sfun=fespro(sig, p, dis, alpha, delta)
        % a = intial estimate of sigma divided by 100
        z = dis/sig;
        sfun = mean(rhork2(z, p, alpha)) - delta;
        
    end

    function zz=rhork2(t, p, alpha)
        % Rocke rho function
        xx = t;
        z  = chi2inv(1-alpha, p);
        % g is \gamma in equation (6.40) of MMY
        g = min(z/p - 1, 1);
        uu1 = xx <= 1-g;
        uu2 = xx > 1-g & xx <= 1+g;
        uu3 = xx > 1+g;
        zz = xx;
        x1 = xx(uu2);
        if ~isempty(x1)
            dd = ((x1-1)/(4*g)).*(3-((x1-1)/g).^2)+.5;
            zz(uu2) = dd;
        end
        
        zz(uu1) = 0;
        zz(uu3) =1;
    end

    function  zz=wrk2(t, p, alpha)
        % Rocke w function (this is the derivative of Rocke rho function
        % This is 
        xx = t;
        z  = chi2inv(1-alpha, p);
        g = min(z/p - 1, 1);
        uu1 = xx <= (1-g);
        uu2 = xx > 1-g & xx <= 1+g;
        uu3 = xx > (1+g);
        zz = xx;
        x1 = xx(uu2);
        dd = (3/(4*g))*(1-((x1-1)/g).^2);
        zz(uu1)=0;
        zz(uu2)=dd;
        zz(uu3)=0;
    end


    function outdes=descen1(x, mu, v, alpha, sig0, dis0)
        
        % one step of reweighted
        [n,p] = size(x);
        
        p2 =-(1/p);
        delta = 0.5*(1-p/n);
        daux= dis0/sig0;
        w = wrk2(daux, p, alpha);
        
        % newloc = new estimate of location using the weights previously found
        % newloc = \sum_{i=1}^n y_i w(d_i) / \sum_{i=1}^n w(d_i)
        mui = sum(bsxfun(@times,x,w),1)/sum(w);
        % Res = n x v matrix which contains deviations from the robust estimate
        % of location
        Res = bsxfun(@minus,x, mui);
        % va = new shape matrix
        va= (Res')*bsxfun(@times,Res,w);
        
        disires = disesca(x,mui,va,alpha,delta);
        disi = disires.dis;
        sigi = disires.sig;
        deti = disires.det;
        vi =  disires.v;
        outdes.mui=mui;
        outdes.vi=vi;
        outdes.sig=sigi;
        outdes.dis=disi;
    end

   function outgra=gradien(x, mu1, v1, mu0, v0, sig0, dis0, alpha)
          
        [n,p]=size(x);
        delta = 0.5*(1-p/n);
        % decreG = 0;
        muminG =mu0;
        vminG = v0;
        disminG = dis0;
        sigminG = sig0;
        for i=1:10
            gamma = i/10;
            mu = (1-gamma)*mu0+gamma*mu1;
            v = (1-gamma)*v0+gamma*v1;
            disresf = disesca(x,mu,v,alpha,delta);
                sigma = disresf.sig;
                % decreG = (sig0-sigma)/sig0;
                if (sigma < sigminG) 
                    muminG = mu;
                    vminG = disresf.v;
                    sigminG = disresf.sig;
                    disminG = disresf.dis;
                end
        end
       outgra=struct;
       outgra.mu=muminG;
       outgra.v=vminG;
       outgra.sig=sigminG;
       outgra.dis=disminG;
    end

end