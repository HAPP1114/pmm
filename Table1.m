clc; clear; close all; warning off
addpath(genpath(fileparts(mfilename('fullpath'))));

p = 1000;           % number of features 
n = 300;            % number of observation 
K = 10; 		    % number of non zero elemet in  real solution 
ratio = 2;          % range of value in x (10^ratio)
sigma = 0.1;        % noise standard deviation
xekind=1;           % type of real solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opts
pkind=3;%3:MCP,4:SCAD
disp('MCP model is running ...')
opts.pkind = pkind;
opts.tau = 2.7;
opts.N = 100;
opts.Lmax  = 1;
opts.Lmin = 1e-10;
opts.maxit =5; 
opts.mu = n/log(p);
opts.reltol = 1e-6;
opts.beta0 = zeros(p,1);
opts.sel = 'hbic';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All available solvers:
FIELDS={'PMM'};
testnum=100;
nmethod=1;
CC=0.2;
len=length(CC);
Time =  zeros(testnum,1); 
RE =  zeros(testnum,1);
CM=  zeros(testnum,1);
MS=  zeros(testnum,1);
seednum=0;
beltaes=zeros(testnum,p);
for k=1:len
    cor=CC(k);
    for ii=1:testnum
        seednum = seednum+k*ii;        % the seed number 
        [X,Xt,D,b,be,xe,Ae] = gendata(n,p,K,sigma,ratio,seednum,xekind,cor);
        opts.del = norm(b-be);
        opts.xe=xe;
        opts.D=D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PMM
        ff='PMM';
        if ismember(ff,FIELDS)  
           disp('PMM is runing....')
           tic,
           [beltapmm,lamm,outputpmm] = PMM(X,Xt,b,opts,opts.pkind);
           pmmtime = toc;
           Time(ii,k,1)=pmmtime;
           beltapmm = D*beltapmm;
           beltaes(ii,:)=beltapmm';
           RE(ii)=norm(beltapmm-xe)/norm(xe);
           Apmm=find(beltapmm);
           MS(ii)=length(Apmm);
           if MS(ii)==length(Ae)
              if Apmm==Ae
                 CM(ii)=1; 
              end
           end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output
meanbelta=mean(beltaes);
meanbelta1=meanbelta(Ae);
stdbelta=std(beltaes);
stdbelta1=stdbelta(Ae);
%RE
meanre=mean(RE);
stdre=std(RE);
%MS
meanms=mean(MS);
stdms=std(MS);
%CM
meancm=mean(CM);
stdcm=std(CM);





