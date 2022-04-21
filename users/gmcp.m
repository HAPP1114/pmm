function y = gmcp(beta1,lambda,tau)
p=length(beta1);
y=zeros(p,1);
for i=1:p
    if abs(beta1(i))>lambda*tau
        y(i)=lambda*sign(beta1(i));
    else
        y(i)=beta1(i)/tau;
    end
end
end

