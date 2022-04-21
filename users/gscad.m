function y = gscad(beta1,lambda,tau)
p=length(beta1);
y=zeros(p,1);
for i=1:p
    if abs(beta1(i))>lambda*tau
        y(i)=lambda*sign(beta1(i));
    elseif abs(beta1(i))<=lambda
        y(i)=0;
    else
        y(i)=(sign(beta1(i))*(abs(beta1(i))-lambda))/(tau-1);
    end
end
end

