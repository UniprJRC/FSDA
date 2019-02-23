function [smo,span]=rlsmo(x,y,w,span)
% rlsmo computes a running-lines smoother with global cross-validation

% note that x is sorted

n=length(x);
cross=0;
if span == 0
    cross=1;
end
penal=0.01;
cvmin=1e15;

if cross==1
    cvspan=[0.3,0.4,0.5,0.6,0.7,1.0];

    % cvrss = Matrix which will contain the cross validation residual sum
    % of squares for each value of vector cvspan
    cvrss=zeros(6,1);
    % SMO = matrix which contains smoothed values for each of the 6 values of
    % cvspan
    SMO=zeros(n,6);
    
    for  k=1:6
        % [smo,rss]=smth(x,y,w,span,dof,n,cross);
        disp(['k=' num2str(k)])
        [smo,cvrss(k)]=smth(x,y,w,cvspan(k), n,1);
        
        % smoCHK=smooth(x,y,'lowess',cvspan(k))
        SMO(:,k)=smo;
        if cvrss(k)>cvmin
            % condwhile=0;
        else
            cvmin=cvrss(k);
            idmin=k;
        end
        %subroutine rlsmo
    end
    
    % span=cvspan(idmin);
    if penal <=0
        return
    else
        cvmin= (1.+penal)*cvmin;
        
        jk=7;
        for k=6:-1:1
            jk=jk-1;
            if cvrss(k) < cvmin
                break
            end
        end
        
        span=cvspan(jk);
    end
% else
%     disp('ciao')
end

% Call smth function with cross validation parameter set to 0
[smo,rss,dof,meay]=smth(x,y,w,span,0);
smo=smo+meay;
end


function [smo,rss,dof,meay] = smth(x,y,w,span,cross)
% smoothing function for aperiodic data, uses weights.

% if nargin<5
%     vsmlsq=0;
% end

dof=0;
s0=0;
% Get the dimensions of the input.
n = length(y);

wy = w.*y;
wxy = wy.*x;
data = [cumprod([w,x,x],2),wy,wxy];

% FS: The alternative code is
% data =[w,x, x.^2, wy, wxy]
sw=sum(w);
meay=sum(y.*w)/sw;


if span <1
    ispan=n.*span;
    m=fix(ispan./2);
    
    if m<1
        m=1;
    end
    k = 2*m + 1;
    
    % Compute sum(w), sum(w*x), sum(w*x^2), sum(w*y), sum(w*x*y) over k points.
    sums = zeros(n,5);
    % ----------------------------------------------
    % Slower, more accurate code:
    % temp = filter(ones(it,1),1,data);
    % sums(m+1:n-m,:) = temp(k:end,:);
    % ----------------------------------------------
    % Faster, slightly less accurate code:
    cs = [0 0 0 0 0;cumsum(data)];
    sums(m+1:n-m,:) = cs(k+1:end,:) - cs(1:end-k,:);
    % ----------------------------------------------
    % Repeat row m+1 in rows 1:m
    % and repeat row n-m in the last m rows
    sums(1:m,:) = sums((m+1)*ones(1,m),:);
    sums(n-m+1:end,:) = sums((n-m)*ones(1,m),:);
    
    remfromtop=flipud(cumsum(flipud(data(m+2:2*m+1,:)),1));
    sums(1:m,:)=sums(1:m,:)-remfromtop;
    
    remfrombottom=(cumsum((data(n-2*m:n-m-1,:)),1));
    sums(n-m+1:end,:)=sums(n-m+1:end,:)-remfrombottom;
    
else
    sums=repmat(sum(data,1),n,1);
end

% To remove from each row of sums observation i (that is to do cross
% validation)
if cross==1
    sums=sums-data;
end

denom = sums(:,1).*sums(:,3) - sums(:,2).^2;
a = (sums(:,4).*sums(:,3) - sums(:,2).*sums(:,5))./denom;
b = (sums(:,1).*sums(:,5) - sums(:,2).*sums(:,4))./denom;
smo = a + b.*x;
if span<1
    if m==1
        smo(1)=y(2);
        smo(end)=y(end-1);
    end
    
    % Check the presence of NaN inside smo
    % NaN are due to constant values of x over the span
    NaNsmo=isnan(smo);
    
    % The smoothed values are simply equal to the weighted average of y over
    % the span, if x is constant over the span
    if sum(NaNsmo)>0
        smo(NaNsmo)=sums(NaNsmo,4)./sums(NaNsmo,1);
    end
    % Return smoothed values in terms of deviation from the overall mean of y
    smo=smo-meay;
else
    smo=smo-sums(:,4)./sums(:,1);
end

rss=sum((y-meay-smo).^2.*w)/sw;

% rss=0;
%
% rss=0.0;
% for  i=1:n
%     rss=rss+(w(i)./sumw).*(y(i)-s0-smo(i)).^2;
% end

end



