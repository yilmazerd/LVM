function [Err,Jake]=LVcompute(x,Tmeas,row,col,scale,fig1,fig2);

x=x.';
Tmod=zeros(row,col);

%For Gradient based
C=x(1:row^2);
Cmat=zeros(row);
Cmat(:)=C/scale;
Cmat=Cmat.';
A=x(row^2+1:length(x));
EA=diag(exp(A),0);
INVA=diag(1./A,0);
EAm1=EA-eye(row);
Tmod(:,1)=Tmeas(:,1);
for ind=2:col
    T0=Tmod(:,ind-1);
    num=EA*T0;
    den=ones(row,1)-INVA*EAm1*Cmat*T0;
    Tmod(:,ind)=num./den;
end
ErrM=(Tmod(:,2:col)-Tmeas(:,2:col))./Tmeas(:,2:col);
Err=ErrM(:);

%Tmod

xvars=row+row^2;
dC=zeros(row,row);

JT=zeros(row,xvars);
dTbuffer=zeros(row,xvars);
dparam=zeros(xvars,1);
dparam(1)=1/scale;
Jake=zeros(row*(col-1),xvars);

for per=2:col
    T=Tmod(:,per-1);
    dTmod=dTbuffer;
    for ii=1:row
        Tii=T(ii);
        dTii=dTmod(ii,:);
        Aii=A(ii);
        N=exp(Aii)*Tii;
        D=1-(dot(Cmat(ii,:),T)*(exp(Aii)-1)/Aii);
        X=(exp(Aii)-1)/Aii;
        for jj=1:xvars
            dC(:)=dparam(1:row^2);
            dC=dC.';
            dTiijj=dTii(jj);
            dTjj=dTmod(:,jj);
            if(jj==ii+row^2)
                dAi_dp=1;
                dXij_dp=1/Aii*exp(Aii)-(exp(Aii)-1)/(Aii)^2;
            else
                dAi_dp=0;
                dXij_dp=0;
            end 
            dNij_dp=exp(Aii)*dTiijj+Tii*exp(Aii)*dAi_dp;
            dDij_dp=-1*(dXij_dp*dot(Cmat(ii,:),T)+X*(dot(dC(ii,:),T)+dot(Cmat(ii,:),dTjj)));
            JT(ii,jj)=(1/D*dNij_dp-N/D^2*dDij_dp);
            JTnorm(ii,jj)=JT(ii,jj)/Tmeas(ii,per);
            %JTnorm(ii,jj)=JT(ii,jj);
            dparam=circshift(dparam,1);
        end        
    end
    dTbuffer=JT;
    ji=per-1;
    Jake((ji-1)*row+1:ji*row,1:xvars)=JTnorm;
end

