close all
clear all
clc
%Generisanje signala
N=1000;
Fs=100;
n=1:N;
F1=1/100;
F2=3/100;
F3=7/100;
x=cos(2*pi*F1*n)+0.5*cos(2*pi*F2*n)+3*cos(2*pi*F3*n);
h=power_line_noise_filter(100,7,30,0.5);
b=h(:,1);
a=h(2:end,2);
b=b';
a1=h(:,2);
a1=a1';
a=a';
%filtriranje signala i crtanje spektara
y=IIR_direct_II(b,a,x);
nfft=1024*4;
X=abs(fft(x,nfft));
f=(0:1/nfft:((nfft-1)/nfft))*Fs;
Y=abs(fft(y,nfft));
figure
plot(f,X),xlim([0 Fs/2]),title('Spektar signala x[n]'),ylabel('|X(jf)|'),xlabel('Frekvencija [Hz]'), grid on;
figure
plot(f,Y),xlim([0 Fs/2]),title('Spektar signala y[n]'),ylabel('|Y(jf)|'),xlabel('Frekvencija [Hz]'), grid on;
%impulsni odziv filtra i njegovo crtanje
figure
impuls=[1  zeros(1,999)];
h1=filter(b,a1,impuls);
stem(0:49,h1(1:50)),title('Impulsni odziv filtra'), xlabel('n');
%Kada je stabilan filtar idalje
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',58,'ProductFractionLength',41,...
        'SumMode','SpecifyPrecision','SumWordLength',58,'SumFractionLength',41);
fi_params=struct('FILTER_COEFITIENTA_BITLENGTH',16,'FILTER_COEFITIENTA_FRAC',11,...
                   'FILTER_COEFITIENTB_BITLENGTH',21,'FILTER_COEFITIENTB_FRAC',15,...
                   'SIGNAL_BITLENGTH',21,'SIGNAL_FRAC',15);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction='Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
%Odredjivanje nula i polova kad je idalje stabilan filtar
nulefpstab=roots(double(FI_b))
polovifpstab=roots([1 double(FI_a)])
%Isto to samo sad vise nije stabilan
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',57,'ProductFractionLength',40,...
        'SumMode','SpecifyPrecision','SumWordLength',57,'SumFractionLength',40);
fi_params=struct('FILTER_COEFITIENTA_BITLENGTH',15,'FILTER_COEFITIENTA_FRAC',10,...
                   'FILTER_COEFITIENTB_BITLENGTH',21,'FILTER_COEFITIENTB_FRAC',15,...
                   'SIGNAL_BITLENGTH',21,'SIGNAL_FRAC',15);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction = 'Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
%Odredjivanje nula i polova kad vise nije stabilan
nulefpnstab=roots(double(FI_b))
polovifpnstab=roots([1 double(FI_a)])
%Crtanje nula i polova kad je filtar stabilan i kad nije
figure
subplot(122)
[hznstab,hpnstab,htnstab]=zplane(nulefpnstab,polovifpnstab); 
set(findobj(hznstab,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(hpnstab,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(htnstab,'Type','line'),'LineWidth',2);
title('Raspored nula i polova fixed point kad vise nije stabilan')
xlabel('Re(z)');
ylabel('Im(z)');
subplot(121)
[hzstab,hpstab,htstab]=zplane(nulefpstab,polovifpstab); 
set(findobj(hzstab,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(hpstab,'Type','line'),'LineWidth',3,'MarkerSize',10);
set(findobj(htstab,'Type','line'),'LineWidth',2);
title('Raspored nula i polova fixed point kada je idalje stabilan')
xlabel('Re(z)');
ylabel('Im(z)');
%Amplitudska karakteristika koja ne odstupa previse
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',62,'ProductFractionLength',43,...
        'SumMode','SpecifyPrecision','SumWordLength',62,'SumFractionLength',43);
fi_params = struct('FILTER_COEFITIENTA_BITLENGTH',20,'FILTER_COEFITIENTA_FRAC',13,...
                   'FILTER_COEFITIENTB_BITLENGTH',21,'FILTER_COEFITIENTB_FRAC',15,...
                   'SIGNAL_BITLENGTH',21,'SIGNAL_FRAC',15);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction = 'Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
%Kada se filtrira x sa minimalnim brojem bitova za predstavljanje
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',51,'ProductFractionLength',32,...
        'SumMode','SpecifyPrecision','SumWordLength',51,'SumFractionLength',32);
fi_params = struct('FILTER_COEFITIENTA_BITLENGTH',20,'FILTER_COEFITIENTA_FRAC',13,...
                   'FILTER_COEFITIENTB_BITLENGTH',21,'FILTER_COEFITIENTB_FRAC',15,...
                   'SIGNAL_BITLENGTH',10,'SIGNAL_FRAC',4);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction = 'Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(FI_x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
%Filtriranje ulaznog signala
y_fixed_point=FI_IIR_direct_II(FI_b,FI_a,FI_x);
[H,w]=freqz(b,a1,1024*4);
Ha=abs(H);
[Hfp,w]=freqz(double(FI_b),double([1 FI_a]),1024*4);
Hafp=abs(Hfp);
figure
plot(w,Ha,'LineWidth',2),title('Amplitudska karakteristika koja ne odstupa previse'),grid on, hold on,
plot(w,Hafp,'r','LineWidth',2),
xlabel('w');
ylabel('|H(e(jw))|');
%Odredjivanje gresaka i crtanje grafika za razlicite metode
ey=y-double(y_fixed_point);
figure
subplot(411);
plot(n,x);
title('Ulazni signal');
subplot(412);
plot(n,y);
title('Izlazni signal sa double preciznoscu');
subplot(413);
plot(n,y_fixed_point);
title('Izlazni signal sa fixed point preciznoscu');
subplot(414);
plot(n,ey);
title('Razlika izlaznih signala sa double i sa fixed point preciznoscu');
%Prvi primer kad je velika preciznost svih koeficijenata i ulaznog signala
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',135,'ProductFractionLength',98,...
        'SumMode','SpecifyPrecision','SumWordLength',135,'SumFractionLength',98);
fi_params = struct('FILTER_COEFITIENTA_BITLENGTH',55,'FILTER_COEFITIENTA_FRAC',36,...
                   'FILTER_COEFITIENTB_BITLENGTH',55,'FILTER_COEFITIENTB_FRAC',46,...
                   'SIGNAL_BITLENGTH',25,'SIGNAL_FRAC',16);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction = 'Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(FI_x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
%Filtriranje u tom slucaju i crtanje amplitudske karakteristike
y_fixed_point=FI_IIR_direct_II(FI_b,FI_a,FI_x);
[H,w]=freqz(b,a1,1024*4);
Ha=abs(H);
[Hfp,w]=freqz(double(FI_b),double([1 FI_a]),1024*4);
Hafp=abs(Hfp);
figure
plot(w,Ha,'LineWidth',2),title('Amplitudska karakteristika kad je velika preciznost'),grid on, hold on,
plot(w,Hafp,'r','LineWidth',2),
xlabel('w');
ylabel('|H(e(jw))|');
%Odredjivanje greske i crtanje
ey=y-double(y_fixed_point);
figure
subplot(411);
plot(n,x);
title('Ulazni signal');
subplot(412);
plot(n,y);
title('Izlazni signal sa double preciznoscu');
subplot(413);
plot(n,y_fixed_point);
title('Izlazni signal sa fixed point preciznoscu');
subplot(414);
plot(n,ey);
title('Razlika izlaznih signala sa double i sa fixed point preciznoscu');
%Drugi primer kad je manja preciznost koeficijenata u imeniocu
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',105,'ProductFractionLength',78,...
        'SumMode','SpecifyPrecision','SumWordLength',105,'SumFractionLength',78);
fi_params = struct('FILTER_COEFITIENTA_BITLENGTH',25,'FILTER_COEFITIENTA_FRAC',16,...
                   'FILTER_COEFITIENTB_BITLENGTH',55,'FILTER_COEFITIENTB_FRAC',46,...
                   'SIGNAL_BITLENGTH',25,'SIGNAL_FRAC',16);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction = 'Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(FI_x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
%Filtriranje u tom slucaju
y_fixed_point=FI_IIR_direct_II(FI_b,FI_a,FI_x);
[H,w]=freqz(b,a1,1024*4);
Ha=abs(H);
[Hfp,w]=freqz(double(FI_b),double([1 FI_a]),1024*4);
Hafp=abs(Hfp);
%Crtanje amplitudske karakteristike u tom slucaju
figure
plot(w,Ha,'LineWidth',2),title('Amplitudska karakteristika kad je manja preciznost koeficijenata u imeniocu'),grid on, hold on,
plot(w,Hafp,'r','LineWidth',2),
xlabel('w');
ylabel('|H(e(jw))|');
%Odredjivanje i crtanje greske
ey=y-double(y_fixed_point);
figure
subplot(411);
plot(n,x);
title('Ulazni signal');
subplot(412);
plot(n,y);
title('Izlazni signal sa double preciznoscu');
subplot(413);
plot(n,y_fixed_point);
title('Izlazni signal sa fixed point preciznoscu');
subplot(414);
plot(n,ey);
title('Razlika izlaznih signala sa double i sa fixed point preciznoscu');
%Treci primer kad je manja preciznost koeficijenata u brojiocu
FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Saturate',...
        'ProductMode','SpecifyPrecision','ProductWordLength',105,'ProductFractionLength',78,...
        'SumMode','SpecifyPrecision','SumWordLength',105,'SumFractionLength',78);
fi_params = struct('FILTER_COEFITIENTA_BITLENGTH',55,'FILTER_COEFITIENTA_FRAC',46,...
                   'FILTER_COEFITIENTB_BITLENGTH',25,'FILTER_COEFITIENTB_FRAC',16,...
                   'SIGNAL_BITLENGTH',25,'SIGNAL_FRAC',16);
FI_b=fi(b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
FixedPointAttributes.OverflowAction = 'Wrap';
FI_b=fi(FI_b,true,fi_params.FILTER_COEFITIENTB_BITLENGTH,fi_params.FILTER_COEFITIENTB_FRAC,FixedPointAttributes);
FI_a=fi(FI_a,true,fi_params.FILTER_COEFITIENTA_BITLENGTH,fi_params.FILTER_COEFITIENTA_FRAC,FixedPointAttributes);
FI_x=fi(FI_x,true,fi_params.SIGNAL_BITLENGTH,fi_params.SIGNAL_FRAC,FixedPointAttributes);
%Filtriranje u tom slucaju i crtanje amplitudske karakteristike
y_fixed_point=FI_IIR_direct_II(FI_b,FI_a,FI_x);
[H,w]=freqz(b,a1,1024*4);
Ha=abs(H);
[Hfp,w]=freqz(double(FI_b),double([1 FI_a]),1024*4);
Hafp=abs(Hfp);
figure
plot(w,Ha,'LineWidth',2),title('Amplitudska karakteristika kad je manja preciznost koeficijenta u brojiocu'),grid on, hold on,
plot(w,Hafp,'r','LineWidth',2),
xlabel('w');
ylabel('|H(e(jw))|');
%Odredjivanje greske i crtanje 
ey=y-double(y_fixed_point);
figure
subplot(411);
plot(n,x);
title('Ulazni signal');
subplot(412);
plot(n,y);
title('Izlazni signal sa double preciznoscu');
subplot(413);
plot(n,y_fixed_point);
title('Izlazni signal sa fixed point preciznoscu');
subplot(414);
plot(n,ey);
title('Razlika izlaznih signala sa double i sa fixed point preciznoscu');




