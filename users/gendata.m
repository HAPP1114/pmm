function [X,Xt,D,y,ye,xe,Ae]=gendata(n,p,K,sigma,ratio,seednum,xekind,cor)
%========================================================================
% INPUTS:                                                                
%     n   ---- number of samples                                         
%     p   ---- signal length                                            
%     K   ---- number of nonzero elements in the signal                  
%  ratio  ---- range of value in the signal (= 10^ratio)                 
%  sigma  ---- noise variance                                                                              
%    cor  ---- corr. coeff.                              
% OUTPUTS:                                                               
%     X   ---- normalized sample matrix     
%     Xt  ---- transpose 
%     D   ---- diagnal mtrix makes X normalized                          
%     y   ---- data vector with noise                                    
%    ye   ---- data vector without noise                                 
%    xe   ---- true signal                                               
%     Ae   ---- the support of xe                                                                                   %
%========================================================================

rand('seed',seednum);   % fix seed
randn('seed', seednum); % fix seed

% generate signal
xe = zeros(p,1);     % true signal  
switch xekind
    case 1
        %fixed case for Table 1
        Ae=[30;198;269;395;442;495;637;766;777;865];
        xe(Ae)=[6;-11;-10;25;-8;100;-9;-10;5;1];
    case 2
        %random case
        q = randperm(p);
        Ae = q(1:K);
        if ratio ~= 0
           vx = ratio*rand(K,1);
           vx = vx-min(vx);
           vx = vx/max(vx)*ratio;
           xe(Ae) = 10.^vx.*sign(randn(K,1));
        else
           xe(Ae) = sign(randn(K,1));
        end
end

% generate matrix X
X = randn(n,p);
if cor ~= 0
   Sigma = zeros(p,p);
   for k = 1:p
       for l = 1:p 
           Sigma(k,l)=cor^(abs(k-l));
       end
   end
   X= X*chol(Sigma);
end

% generate right hand side
ye  = X*xe;
y   = ye + sigma*randn(n,1);
[X,D] = normalize(X); 
Xt=X';
end