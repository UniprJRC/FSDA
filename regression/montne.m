function ymon=montne(y)
% montne forces input vector y to be monotonic
%
%   This routine returns a vector of values that minimizes the
%   sum of squares \sum_{i=1}^n (y - \hat y ).^2 under the monotonicity constraint that
%   $\hat y(i) \geq \hat y(j)$ for $i>j$, i.e., the values in $\hat y$ are
%   monotonically non-decreasing. This routine uses the "pool adjacent
%   violators" algorithm (PAVA).



%% Beginning of code

% Initialize fitted values (monotonic y) to the given values.
yhat = y;
n = length(yhat);
w=ones(1,n);

% Written by Maigo on 8/14/2012 to reduce the complexity from O(n^2) to O(n)
b = 0; 
bstart = zeros(1,n); 
bend = bstart;
for i = 1:n
    b = b + 1;
    yhat(b) = yhat(i);
    w(b) = w(i);
    bstart(b) = i; bend(b) = i;
    while b > 1 && yhat(b) < yhat(b-1)
        yhat(b-1) = (yhat(b-1) * w(b-1) + yhat(b) * w(b)) / (w(b-1) + w(b));
        w(b-1) = w(b-1) + w(b);
        bend(b-1) = bend(b);
        b = b - 1;
    end
end
idx = zeros(1,n);
for i = 1:b
    idx(bstart(i) : bend(i)) = i;
end
% Maigo end

% ymon = Monotonic (non decreasing)  fitted values  
ymon=yhat(idx);
end

% The original routine literarly translated from pure Fortran 77 is given below
% 
% function xmon=montne(x,n)
% % Apply isothonic regression to vector x
% bb=0;
% eb=bb;
% 
% while eb<=n-1
%     bb=eb+1;
%     eb=bb;
%     if eb<n-1
%         
%         vvv=1;
%         while vvv==1
%             if (x(bb) == x(eb+1))
%                 if eb==n-1
%                     vvv=0;
%                 else
%                     eb=eb+1;
%                 end
%             else
%                 vvv=0;
%             end
%         end
%     else
%     end
%     
%     if eb>=n
%         jjj=0;
%     else
%         jjj=1;
%     end
%     
%     while jjj==1
%         if eb<=n-1 && (x(eb) > x(eb+1)) % monotonicity is violated
%             br=eb+1;
%             er=br;
%             
%             zzz=1;
%             while zzz==1
%                 if er>=n
%                     pmn=(x(bb)*(eb-bb+1)+x(br)*(er-br+1))/(er-bb+1);
%                     zzz=0;
%                 else
%                     if (x(er+1) == x(br))
%                         er=er+1;
%                     else
%                         pmn=(x(bb)*(eb-bb+1)+x(br)*(er-br+1))/(er-bb+1); % take average
%                         zzz=0;
%                     end
%                 end
%             end
%             
%             eb=er;
%             x(bb:eb)=pmn;
%         else
%         end
%         
%         
%         if bb<=1
%             break
%         else
%             bbless1=0;
%         end
%         
%         
%         if bbless1==1
%             continue
%         end
%         
%         if    x(bb-1) > x(bb)
%             bl=bb-1;
%             el=bl;
%             while bl>1
%                 if (x(bl-1)~= x(el))
%                     break   
%                 else
%                     bl=bl-1;
%                 end
%             end
%             pmn=(x(bb)*(eb-bb+1)+x(bl)*(el-bl+1))/(eb-bl+1);
%             bb=bl;
%             x(bb:eb)=pmn;
%         else
%             break
%         end
%     end
%     
% end
% 
% % Last iteration 
% 
% if bb>1 && x(bb-1) > x(bb)
%     bl=bb-1;
%     el=bl;
%     while bl>1
%         if (x(bl-1)~= x(el))
%             break
%         else
%             bl=bl-1;
%         end
%     end
%     pmn=(x(bb)*(eb-bb+1)+x(bl)*(el-bl+1))/(eb-bl+1);
%     bb=bl;
%       x(bb:eb)=pmn;
%       
%     while bb>1
%         if    x(bb-1) > x(bb)
%             bl=bb-1;
%             el=bl;
%             while bl>1
%                 if (x(bl-1)~= x(el))
%                     break  
%                 else
%                     bl=bl-1;
%                 end
%             end
%             pmn=(x(bb)*(eb-bb+1)+x(bl)*(el-bl+1))/(eb-bl+1);
%             bb=bl;
%               x(bb:eb)=pmn;
%         else
%             if eb>=n
%                 break
%             end
%         end
%     end
% end
% 
% xmon=x;
% end
% 
