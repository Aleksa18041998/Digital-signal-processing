close all
clear all
clc
%Ucitavanje signala i plotovanje
ecg=load('ecg_corrupted.mat');
ecg1=struct2array(ecg);
fs=360;
Ts=1/fs;
t=0:Ts:((length(ecg1)-1))*Ts;
figure
plot(t,ecg1),grid on;
title('EKG signal');
xlabel('Vreme [s]');
%Filtriranje VF filtrom i plotovanje isfiltriranog signala
fa1=0.4;
fp1=1;
Aa1=30;
Ap1=0.5;
y1=baseline_drift_filter(fs,fa1,fp1,Aa1,Ap1);
b1=y1(:,1);
a1=y1(:,2);
ecgbasefiltered=filter(b1,a1,ecg1);
t=0:Ts:((length(ecgbasefiltered)-1))*Ts;
figure
plot(t,ecgbasefiltered),grid on;
title('EKG signal filtriran VF');
xlabel('Vreme [s]');
%Filtriranje NO filtrom i plotovanje konacnog signala
fc2=60;
Aa2=40;
Ap2=0.5;
y2=power_line_noise_filter(fs,fc2,Aa2,Ap2);
b2=y2(:,1);
a2=y2(:,2);
ecgpowerfiltered=filter(b2,a2,ecgbasefiltered);
t=0:Ts:((length(ecgpowerfiltered)-1))*Ts;
figure
plot(t,ecgpowerfiltered),grid on;
title('EKG signal filtriran NO');
xlabel('Vreme [s]');