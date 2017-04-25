function Tout=LVTime(s,c,a,T0,years);

Tbuffer=zeros(length(a),years+1);
Tbuffer(:,1)=T0;

for ind=1:length(Tbuffer)
    T=Tbuffer(:,ind);
    scprod=s.*c;
    scale=(exp(a)-1)./a;
    den=1-(scprod*T).*scale;
    num=exp(a).*T;
    Tbuffer(:,ind+1)=num./den;
end

Tout=Tbuffer;
disp(Tout');