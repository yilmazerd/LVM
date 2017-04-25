clear all;close all;
% pwood = [1	1.2	1.3	1.34 1.4 1.5 1.6];
% wbo = [0.14	0.24	0.3308	0.3327 0.4 0.5 0.7];
% wbo2 = [0.02 0.05 0.1 0.12 0.14 0.15 0.17];
% pwood = [18.0617	15.8378	17.1886	18.539	19.4882	19.72	19.6069	18.2317	16.2252	17.0305	19.5575	19.6745	20.8815	22.0313	21.6602	21.0311	19.541	17.2763	17.7375	17.6534];
% wbo = [0.14	0.24	0.3308	0.3327	0.349	0.4513	0.6256	0.7713	0.9883	1.58	2.42	2.8451	3.68	4.6857	5.2759	6.3396	6.8007	6.6162	8.34	9.4606];
% wbo2 = wbo + wbo*19*rand();

BE =[9074 12795 46832 60368 64175];
PE =[6966 37558 41376 56548 49118];
% BEPE = BE/2+PE/2;
% BE=BEPE;
% PE=BEPE;
HV =[269178 434498 495685 452152 384417];
firstyear = 2011;
lastyear = 2017;

% BE = [0.1 0.15 0.25 0.35 0.65 1 1.7 3 7 16]
% HV = [8.0 8.2 8.3 8.4 8.5 8.6 8.7 8.8 8.82 8.83]
% Data = [BE;PE;HV];
Data = [HV;BE;PE];
a = rand(3,1);
% s = [1 1 1;1 1 1; 1 1 1]'; %important to have transpose
s = ones(3,3) - rand(3,3)*2;
%s = [1 -1;-1 1];
c = [0.00001    0.00001 0.00001;   0.00001 0.00001 0.0014;0.00001 0.00001 0.0014];
c =rand(3,3)/10;

[sout,cout,aout]=LVopt(Data,true,a,s,c);

T0=Data(:,1);
Tout=LVmod(sout,cout,aout,T0,lastyear-firstyear-1);
yr=firstyear:1:lastyear;
figure
plot(yr,Tout,'.-','LineWidth',2);
hold
plot(yr(1:max(size(Data))), Data, 'sq','LineWidth',3)
grid
xlabel('Year')
ylabel('Units')


%% COMPUTE SUM OF SQUARES
[p k] = size(Data);
disp('Total sum of squares')
sumsq = sum(sum((Data(1:end,2:end) - Tout(1:p,2:k)).^2))/((p-1)*k)

