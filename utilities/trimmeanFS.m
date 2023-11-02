function meantrimmed = trimmeanFS(y,alpha)
%trimmeanFS trims a proportion alpha of observations from both ends
%
% let $\alpha \in [0, \, 0.5)$ and $m=[(n-1) \alpha]$ where $[\cdot]$ is the symbol which denotes the integer part,
%the $\alpha-$trimmed mean $(\overline y_\alpha)$ is defined as
% \[
% \overline y_\alpha = \frac{1}{n-2m} \sum_{i=m+1}^{n-m} y_{[i]}.
% \]
n=length(y);
ysor=sort(y);
m=floor((n-1)*alpha);
meantrimmed=mean(ysor(m+1:n-m));
end


