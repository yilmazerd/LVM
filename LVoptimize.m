function [s,c,a]=LVoptimize(T,LIM,a,s,c,sc_Ue,sc_Le)
%[sout,cout,aout]=LVopt(T,LIM,a,s,c)
%Erdem Yilmaz - Chris Aden
%T=T[m,n] is a technology matrix containing m technology rows specified
%LIM: boolean: true-use simple constraints on parameters, false-unconstrained
%minimization
%over n time periods
%s=initialized correlation matrix
%c=initialized dependencies matrix
%a=initialized growth vector
%Uawe ma

[m,n]=size(T);
scale=1e10; %optimization scalar
sc_U=ones(m^2,1);
sc_L=-1*ones(m^2,1);
if nargin==2;
    a_init=ones(m,1);
    sc_init=zeros(m);
elseif nargin==3 && length(a)==m;
    a_init=a;
    sc_init=zeros(m);
elseif nargin==4 && length(a)==m;
    a_init=a;
    sc_init=s;    
elseif nargin==5 && length(a)==m && isequal(size(s),size(ones(m)));
    a_init=a;
    sc_mat=s.*c;
    sc_mat=sc_mat';
    sc_init=sc_mat(:);
    scale=1/min(abs(sc_init));
    sc_init=sc_init*scale;
elseif nargin==7 && length(a)==m && isequal(size(s),size(ones(m)));
    a_init=a;
    sc_mat=s.*c;
    sc_mat=sc_mat';
    sc_init=sc_mat(:);
    scale=1/min(abs(sc_init));
    sc_init=sc_init*scale;  
    if LIM==true
        sc_U=sc_Ue(:);
        sc_L=sc_Le(:);
    end
end

x_init=[sc_init(:);a_init].';

% fig1=figure;
% set(fig1,'position',[100 600 400 400]);
% meas1=plot(T','-');
% set(meas1,'tag','meas1');
% set(gca,'nextplot','add');
% 
% fig2=figure;
% set(fig2,'position',[100 100 400 400]);
% meas2=semilogy(T','-');
% set(meas2,'tag','meas2');
% set(gca,'nextplot','add');
%set constraints

%set simple optimization bounds
if LIM==true
    %sc_U=ones(m^2,1)*1000
    % sc_U=ones(m^2,1)*1
    %sc_U=sc_U-[1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1]';
    %sc_U = [1 -1 -1 -1 -1 1 -1   -1 -1  -1  1 -1 -1 -1 -1 1]';
    %sc_U=zeros(m^2,1)*1000;
    %sc_L=sc_U*-1;
    %sc_L=sc_U;
    %sc_L = sc_U;
    %sc_L=zeros(m^2,1)*1000;
    %sc_U = [-0.9077    -1.0000    -1.0000    100  -0.7028    -0.3453    -1.0000    100  -1.0000    -1.0000    -1.0000    100 100    100   100    100]';
    %sc_L = [-0.9077    -1.0000    -1.0000    -100 -0.7028    -0.3453    -1.0000    -100 -1.0000    -1.0000    -1.0000    -100 -100    -100   -100    -100]';
   
    a_U=ones(m,1)*100;
    %a_U=[1.4598    0.6087    0.0375 100]';
    a_L=ones(m,1)*0;
    %a_L=[1.4598    0.6087    0.0375 0]';
    U=[sc_U;a_U].';
    L=[sc_L;a_L].';
elseif LIM==false
    % remove constraints
    U=[];
    L=[];
else
    error('Constraint variable is not understood.');
end

tolx=1e-8;
tolf=1e-8;
iter=100;
funeval=1e5;


% 2016 Modification
options=optimset('TolX',tolx,'TolFun',tolf,'MaxIter',iter,'Jacobian','on','display','iter',...
    'MaxFunEvals',funeval,'DerivativeCheck','off');

%x=lsqnonlin('LVcomp',x_init,L,U,options,T,m,n,scale,fig1,fig2);
x=lsqnonlin('LVcompute',x_init,L,U,options,T,m,n,scale);


x=x.';
scd=x(1:m^2);
sc=ones(m);
sc(:)=scd;
sc=sc.';
a=x(m^2+1:m^2+m);
s=sign(sc);
c=s.*sc/scale;

a
s
c
x