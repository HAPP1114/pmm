
clc; clear; close all; warning off
addpath(genpath(fileparts(mfilename('fullpath'))));

[label1,instance1]=libsvmread('mpg_scale.txt');
b=label1;
A=instance1;
[m,n] = size(A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Comment out or delete this section when running log1p.E2006
d=7;% order of polynomial
v = mypartition(n+1,d);
v = v';
nn = size(v,2);
AA = zeros(m,nn);
for q=1:nn
    AAq = ones(m,1);
    for j = 1:n
        if v(j,q) > 0
           AAq = AAq.*(A(:,j).^v(j,q));
        end
    end
    AA(:,q) = AAq;
end
A = AA;
n = nn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove zero columns
aa = full(sqrt(sum(A.*A)));
idx = find(aa > 0);
if length(idx) < n
   A = A(:,idx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=A;
Xt=X';
[n,p]=size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vpkind=[3,4];
for kk=1:length(Vpkind)
    pkind=Vpkind(kk);
    switch pkind
       case 3
	     disp('MCP model is running ...')
	     opts.pkind = pkind;
	     opts.tau = 2.7;
	     opts.N = 100;
         opts.Lmax  = 1;
	     opts.Lmin = 1e-10;
	     opts.maxit =5; 
	     opts.mu = n/log(p);
	     opts.reltol = 5e-2;
	     opts.beta0 = zeros(p,1);
	     opts.sel = 'hbic';
       case 4
	     disp('SCAD model is running ...')
	     opts.pkind = pkind;
	     opts.tau = 3.7;
	     opts.N = 100;
         opts.Lmax  = 1;
	     opts.Lmin = 1e-10;
	     opts.maxit = 5; 
	     opts.mu = n/log(p);
	     opts.reltol = 5e-2;
	     opts.beta0 = zeros(p,1);
	     opts.sel ='hbic';		
      otherwise
		error('Undefined penalty !')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All available solvers:
    FIELDS={'PMM','CD2'};
    nmethod=2;
    Time=zeros(nmethod,1);
    NN=zeros(nmethod,1);
    KKT=zeros(nmethod,1);
    MVP=zeros(nmethod,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PMM
    ff='PMM';
    if ismember(ff,FIELDS)  
       disp('PMM is runing....')
       tic,
       [beltapmm,lamm,outputpmm] = PMM(X,Xt,b,opts,opts.pkind);
       pmmtime = toc;
       Time(1)=pmmtime;
       NN(1)=length(find(beltapmm));
       KKT(1)=outputpmm.kkt(end);
       MVP(1)=outputpmm.mvpc1;
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CD
    ff='CD2';
    if ismember(ff,FIELDS) 
       disp('CD is runing....')
       tic,
       [beltacd,lamhat,outputcd]= cdg_path(X,b,opts);
       cdtime = toc;
       Time(2)=cdtime;
       NN(2)=length(find(beltacd));
       KKT(2)=outputcd.KKThat;
       MVP(2)=outputcd.mvpc3;
   end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('NN=');disp(NN);
    disp('KKT=');disp(KKT);
    disp('Time=');disp(Time);
    disp('MVP=');disp(MVP);
    fprintf('tol=%4f, pkind=%4f', opts.reltol, pkind);
end






