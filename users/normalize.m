function [sX,D] = normalize(X)
[n,p] = size(X);
d = 1./sqrt(ones(1,n)*X.^2);
D =  spdiags(d',0,p,p);
sX = X*D;
end