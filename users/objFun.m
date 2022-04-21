function obj = objFun(X,y,beta,lambda,gamma,penalty)
	lost = 0.5*norm(X*beta - y)^2; 
	pen = sum(penFun(beta,lambda,gamma,penalty));
	obj = lost + pen;
end


function pen = penFun(z,lambda,gamma,penalty)
	if penalty == 3
		pen = zeros(size(z));
		idx = abs(z) <= gamma*lambda;
		zz = z(idx);
		pen(idx) = lambda*abs(zz)-abs(zz).^2/(2*gamma);
		idx = abs(z) > gamma*lambda;
		pen(idx) = (1/2)*gamma*(lambda^2);
	else
		pen = zeros(size(z));
		idx = abs(z) <= lambda;
		zz = z(idx);
		pen(idx) = lambda*abs(zz);
		idx = (abs(z) > lambda) & (abs(z) <= gamma*lambda);
		zz = z(idx);
		pen(idx) = -(abs(zz).^2 - 2*gamma*lambda*abs(zz)+lambda^2)/(2*(gamma-1));
		idx = abs(z)>gamma*lambda;
		pen(idx) = (gamma+1)*lambda^2/2;
	end
end
