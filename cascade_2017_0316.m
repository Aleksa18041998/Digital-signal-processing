close all
clear all
clc
%Generisanje ulaznog signala
N=1000;
Fs=100;
n=1:N;
F1=1/100;
F2=3/100;
F3=7/100;
x=cos(2*pi*F1*n)+0.5*cos(2*pi*F2*n)+3*cos(2*pi*F3*n);
h=power_line_noise_filter(100,7,30,0.5);
b=h(:,1);
a=h(:,2);
b=b';
a=a';
%impulsni odziv filtra
figure
impuls=[1  zeros(1,999)];
h1=filter(b,a,impuls);
stem(0:49,h1(1:50)),title('Impulsni odziv filtra'), xlabel('n');
%filtriranje signala i crtanje spektara
sos=tf2sos(b,a);
b=[sos(1,1:3);sos(2,1:3);sos(3,1:3)];
a=[sos(1,5:6);sos(2,5:6);sos(3,5:6)];
y=IIR_direct_II_cascade(b,a,x);
nfft=1024*4;
X=abs(fft(x,nfft));
f=(0:1/nfft:((nfft-1)/nfft))*Fs;
Y=abs(fft(y,nfft));
figure
plot(f,X),xlim([0 Fs/2]),title('Spektar signala x[n]'),ylabel('|X(jf)|'),xlabel('Frekvencija [Hz]'), grid on;
figure
plot(f,Y),xlim([0 Fs/2]),title('Spektar signala y[n]'),ylabel('|Y(jf)|'),xlabel('Frekvencija [Hz]'), grid on;
%Provera da li fixed point funkcija dobro radi i generisanje izlaznog signala i
%njegovo plotovanje
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',57,'ProductFractionLength',45,...
        'SumMode','SpecifyPrecision','SumWordLength',57,'SumFractionLength',45);
fi_params=struct('FILTER_COEFITIENTA_BITLENGTH',19,'FILTER_COEFITIENTA_FRAC',15,...
                 'FILTER_COEFITIENTB_BITLENGTH',19,'FILTER_COEFITIENTB_FRAC',15,...
                 'SIGNAL_BITLENGTH',19,'SIGNAL_FRAC',15);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction = 'Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(FI_x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
y_fixed_point=FI_IIR_direct_II_cascade(FI_b,FI_a,FI_x);
figure
plot(n,y);
title('y cascade');
figure
plot(n,y_fixed_point);
title('y sa fixed point cascade');
%Kad je stabilan idalje filtar koeficijenti u imeniocu
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',19*2+6,'ProductFractionLength',30+3,...
        'SumMode','SpecifyPrecision','SumWordLength',19*2+6,'SumFractionLength',30+3);
fi_params=struct('FILTER_COEFITIENTA_BITLENGTH',6,'FILTER_COEFITIENTA_FRAC',3,...
                   'FILTER_COEFITIENTB_BITLENGTH',19,'FILTER_COEFITIENTB_FRAC',15,...
                   'SIGNAL_BITLENGTH',19,'SIGNAL_FRAC',15);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction='Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
%Nule i polovi kad je idalje filtar stabilan
nulefpstab1=roots(double(FI_b(1,:)))
polovifpstab1=roots([1 double(FI_a(1,:))])
nulefpstab2=roots(double(FI_b(2,:)))
polovifpstab2=roots([1 double(FI_a(2,:))])
nulefpstab3=roots(double(FI_b(3,:)))
polovifpstab3=roots([1 double(FI_a(3,:))])
%Crtanje nula i polova kad je idalje filtar stabilan
figure
subplot(121),
hold on
[hzstab1,hpstab1,htstab1]=zplane(nulefpstab1,polovifpstab1); 
set(findobj(hzstab1,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(hpstab1,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(htstab1,'Type','line'),'LineWidth',2);
[hzstab2,hpstab2,htstab2]=zplane(nulefpstab2,polovifpstab2); 
set(findobj(hzstab2,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(hpstab2,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(htstab2,'Type','line'),'LineWidth',2);
[hzstab3,hpstab3,htstab3]=zplane(nulefpstab3,polovifpstab3); 
set(findobj(hzstab3,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(hpstab3,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(htstab3,'Type','line'),'LineWidth',2);
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
title('Raspored nula i polova fixed point kada je idalje stabilan')
xlabel('Re(z)');
ylabel('Im(z)');
%Isto sad samo kad je postao nestabilan
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',19*2+5,'ProductFractionLength',30+2,...
        'SumMode','SpecifyPrecision','SumWordLength',19*2+5,'SumFractionLength',30+2);
fi_params=struct('FILTER_COEFITIENTA_BITLENGTH',5,'FILTER_COEFITIENTA_FRAC',2,...
                   'FILTER_COEFITIENTB_BITLENGTH',19,'FILTER_COEFITIENTB_FRAC',15,...
                   'SIGNAL_BITLENGTH',19,'SIGNAL_FRAC',15);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction='Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
nulefpnstab1=roots(double(FI_b(1,:)))
polovifpnstab1=roots([1 double(FI_a(1,:))])
nulefpnstab2=roots(double(FI_b(2,:)))
polovifpnstab2=roots([1 double(FI_a(2,:))])
nulefpnstab3=roots(double(FI_b(3,:)))
polovifpnstab3=roots([1 double(FI_a(3,:))])
%Crtanje nula i polova kad je postao nestabilan
subplot(122),
hold on
[hznstab1,hpnstab1,htnstab1]=zplane(nulefpnstab1,polovifpnstab1); 
set(findobj(hznstab1,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(hpnstab1,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(htnstab1,'Type','line'),'LineWidth',2);
[hznstab2,hpnstab2,htnstab2]=zplane(nulefpnstab2,polovifpnstab2); 
set(findobj(hznstab2,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(hpnstab2,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(htnstab2,'Type','line'),'LineWidth',2);
[hznstab3,hpnstab3,htnstab3]=zplane(nulefpnstab3,polovifpnstab3); 
set(findobj(hznstab3,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(hpnstab3,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(htnstab3,'Type','line'),'LineWidth',2);
xlim([-1.6 1.6]);
ylim([-1.6 1.6]);
title('Raspored nula i polova fixed point kada je nestabilan')
xlabel('Re(z)');
ylabel('Im(z)');
%Filtar koji ne odstupa previse
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',43,'ProductFractionLength',31,...
        'SumMode','SpecifyPrecision','SumWordLength',43,'SumFractionLength',31);
fi_params=struct('FILTER_COEFITIENTA_BITLENGTH',12,'FILTER_COEFITIENTA_FRAC',8,...
                   'FILTER_COEFITIENTB_BITLENGTH',12,'FILTER_COEFITIENTB_FRAC',8,...
                   'SIGNAL_BITLENGTH',19,'SIGNAL_FRAC',15);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction='Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
%Racunanje i crtanje amplitudske karakteristike filtra koji ne odstupa
%previse
[H1,w]=freqz(b(1,:),[1 a(1,:)],1024*4);
[H2,w]=freqz(b(2,:),[1 a(2,:)],1024*4);
[H3,w]=freqz(b(3,:),[1 a(3,:)],1024*4);
H=H1.*H2.*H3;
Ha=abs(H);
[H1fp,w]=freqz(double(FI_b(1,:)),double([1 FI_a(1,:)]),1024*4);
[H2fp,w]=freqz(double(FI_b(2,:)),double([1 FI_a(2,:)]),1024*4);
[H3fp,w]=freqz(double(FI_b(3,:)),double([1 FI_a(3,:)]),1024*4);
Hfp=H1fp.*H2fp.*H3fp;
Hafp=abs(Hfp);
figure
plot(w,Ha,'LineWidth',2),title('Amplitudska karakteristika kad fixed point ne odstupa mnogo'),grid on, hold on,
plot(w,Hafp,'r','LineWidth',2),
xlabel('w');
ylabel('|H(e(jw))|');
