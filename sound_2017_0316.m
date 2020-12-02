close all
clear all
clc
%Ucitavanje zagadjenog zvuka i formiranje njegovog spektrograma
[corrupt,Fs]=audioread('sound_corrupted.wav');
window=hamming(512);
nooverlap=256;
[s,f,t]=spectrogram(corrupt,window,nooverlap,256*32,Fs);
S=20*log(abs(s));
figure
imagesc(t,f,S);
title('Corrupted sound');
xlabel('Vreme [s]');
ylabel('Frekvencija [Hz]');
xticks(0:0.2:3);
yticks(0:500:24000);
%Zadavanje parametara za filtar NO1 i odredjivanje prelazne zone 
fp11=1100;
fa11=1300;
fa21=1700;
fp21=1900;
Fp11=fp11/Fs*2;
Fa11=fa11/Fs*2;
Fa21=fa21/Fs*2;
Fp21=fp21/Fs*2;
Bt11=(fa11-fp11)/Fs*2*pi;
Bt21=(fp21-fa21)/Fs*2*pi;
Bt1=min(Bt11,Bt21);
Ap=0.5;
Aa=40;
%Odredjivanje reda filtra
dp=(10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);	
da=10^(-0.05*Aa);
D=(0.005309*log10(dp)*log10(dp)+0.07114*log10(dp)-0.4761)*log10(da);
D=D-(0.00266*log10(dp)*log10(dp)+0.5941*log10(dp)+0.4278);
f=11.01217+0.51244*(log10(dp)-log10(da));
M1=2*pi*D/Bt1-f*Bt1/(2*pi)+1;
M1=ceil(M1);
N1=M1-1;
disp('Pocetni red filtra NO1 je:');
disp(N1);
%projektovanje NO1 filtra optimizacionom metodom
Hd1=[1 1 0 0 1 1];
F1=[0 Fp11 Fa11 Fa21 Fp21 1];
h1=firpm(N1,F1,Hd1);
[H1,w1]=freqz(h1,1,256*32);
Ha1=abs(H1);
figure
%Cranje amplitudske karakteristike filtra koji mozda ne zadovoljava zadate
%gabarite
plot(w1,Ha1),title('Amplitudska karakteristika filtra NO1(firmp)'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
%Crtanje gabarita koje treba da zadovolji
hold on
x1h=[0 Fp11*pi]; y1h=[1-dp 1-dp]; 
x2h=[0 Fp11*pi]; y2h=[1+dp 1+dp];
x1v=[Fp11*pi Fp11*pi]; y1v=[1+dp 0];
x3h=[Fa11*pi Fa21*pi]; y3h=[da da];
x3v=[Fa11*pi Fa11*pi]; y3v=[da 1+dp];
x4v=[Fa21*pi Fa21*pi]; y4v=[da 1+dp];
x4h=[Fp21*pi pi]; y4h=[1+dp 1+dp];
x5h=[Fp21*pi pi]; y5h=[1-dp 1-dp];
x5v=[Fp21*pi Fp21*pi]; y5v=[0 1+dp];
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r',x4v,y4v,'r',x4h,y4h,'r',x5h,y5h,'r',x5v,y5v,'r');
kp11=ceil(256*32*Fp11)+1;
ka11=floor(256*32*Fa11)+1;
ka21=ceil(256*32*Fa21)+1;
kp21=floor(256*32*Fp21)+1;
Hn=Ha1(ka11:ka21)';
Hp=[Ha1(1:kp11)' Ha1(kp21:256*32)'];
%Petlja kojom se proverava da li su zadovoljeni gabariti 
if((max(Hn)<=da)&&(min(Hp)>=(1-dp))&&(max(Hp)<=(1+dp)))
    
else
    while(1)
        N1=N1+1;
        M1=M1+1;
        h1=firpm(N1,F1,Hd1);
        [H1,w1]=freqz(h1,1,256*32); 
        Ha1=abs(H1);
        Hn=Ha1(ka11:ka21)';
        Hp=[Ha1(1:kp11)' Ha1(kp21:256*32)'];
        if((max(Hn)<=da)&&(min(Hp)>=(1-dp))&&(max(Hp)<=(1+dp)));
            break;
        end
    end
end
disp('Konacni red filtra NO1 je:');
disp(N1);
%Crtanje amplitudske karakteristike konacnog NO1 filtra dobijenog
%optimizacionom metodom
figure
plot(w1,Ha1),title('Amplitudska karakteristika ispravljenog filtra NO1(firpm)'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
hold on
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r',x4v,y4v,'r',x4h,y4h,'r',x5h,y5h,'r',x5v,y5v,'r');
%Ciscenje jednog suma u zvuku
y1=conv(corrupt',h1);
%Crtanje spektrograma nakon ciscenja filtrom NO1 
[s1,f1,t1]=spectrogram(y1,window,nooverlap,256*32,Fs);
S1=20*log(abs(s1));
figure
imagesc(t1,f1,S1);
title('Sound ociscen NO1 optimizacionim metodom');
xlabel('Vreme [s]');
ylabel('Frekvencija [Hz]');
xticks(0:0.2:3);
yticks(0:500:24000);
%Zadavanje parametara i gabarita NO2 filtra i odredjivanje prelazne zone
fp12=2200;
fa12=2400;
fa22=2900;
fp22=3100;
Fp12=fp12/Fs*2;
Fa12=fa12/Fs*2;
Fa22=fa22/Fs*2;
Fp22=fp22/Fs*2;
Bt12=(fa12-fp12)/Fs*2*pi;
Bt22=(fp22-fa22)/Fs*2*pi;
Bt2=min(Bt12,Bt22);
Ap=0.5;
Aa=40;
%Odredjivanje reda filtra potrebnog kako bi se zadovoljili gabaritie
dp=(10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);	
da=10^(-0.05*Aa);
D=(0.005309*log10(dp)*log10(dp)+0.07114*log10(dp)-0.4761)*log10(da);
D=D-(0.00266*log10(dp)*log10(dp)+0.5941*log10(dp)+0.4278);
f=11.01217+0.51244*(log10(dp)-log10(da));
M2=2*pi*D/Bt2-f*Bt2/(2*pi)+1;
M2=ceil(M2);
N2=M2-1;
disp('Pocetni red filtra NO2 je:');
disp(N2);
%Projektovanje NO2 filtra optimizacionim metodom na osnovu prethodno
%odredjenog reda filtra
Hd2=[1 1 0 0 1 1];
F2=[0 Fp12 Fa12 Fa22 Fp22 1];
h2=firpm(N2,F2,Hd2);
%Crtanje amplitudske karakteristike filtra koji mozda ne zadovoljava zadate
%gabarite
[H2,w2]=freqz(h2,1,256*32);
Ha2=abs(H2);
figure
plot(w2,Ha2),title('Amplitudska karakteristika filtra NO2(firpm)'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
%Crtanje gabarita
hold on
x1h=[0 Fp12*pi]; y1h=[1-dp 1-dp]; 
x2h=[0 Fp12*pi]; y2h=[1+dp 1+dp];
x1v=[Fp12*pi Fp12*pi]; y1v=[1+dp 0];
x3h=[Fa12*pi Fa22*pi]; y3h=[da da];
x3v=[Fa12*pi Fa12*pi]; y3v=[da 1+dp];
x4v=[Fa22*pi Fa22*pi]; y4v=[da 1+dp];
x4h=[Fp22*pi pi]; y4h=[1+dp 1+dp];
x5h=[Fp22*pi pi]; y5h=[1-dp 1-dp];
x5v=[Fp22*pi Fp22*pi]; y5v=[0 1+dp];
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r',x4v,y4v,'r',x4h,y4h,'r',x5h,y5h,'r',x5v,y5v,'r');
kp12=ceil(256*32*Fp12)+1;
ka12=floor(256*32*Fa12)+1;
ka22=ceil(256*32*Fa22)+1;
kp22=floor(256*32*Fp22)+1;
Hn=Ha2(ka12:ka22)';
Hp=[Ha2(1:kp12)' Ha2(kp22:256*32)'];
%Provera da li su zadovoljeni gabariti i povecavanje reda sve dok se ne
%zadovolje
if((max(Hn)<=da)&&(min(Hp)>=(1-dp))&&(max(Hp)<=(1+dp)))
    
else
    while(1)
        N2=N2+1;
        M2=M2+1;
        h2=firpm(N2,F2,Hd2);
        [H2,w2]=freqz(h2,1,256*32); 
        Ha2=abs(H2);
        Hn=Ha2(ka12:ka22)';
        Hp=[Ha2(1:kp12)' Ha2(kp22:256*32)'];
        if((max(Hn)<=da)&&(min(Hp)>=(1-dp))&&(max(Hp)<=(1+dp)))
            break;
        end
    end
end
disp('Konacni red filtra NO2 je:');
disp(N2);
%Konacna amplitudksa karakteristika
figure
plot(w2,Ha2),title('Amplitudska karakteristika ispravljenog filtra NO2(firpm)'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
hold on
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r',x4v,y4v,'r',x4h,y4h,'r',x5h,y5h,'r',x5v,y5v,'r');
%Ciscenje delimicno ociscenog zvuka filtrom NO2 dobijenim optimizacionim
%metodom
y2=conv(y1,h2);
%Spektrogram idalje delimicno ociscenog zvuka
[s2,f2,t2]=spectrogram(y2,window,nooverlap,256*32,Fs);
S2=20*log(abs(s2));
figure
imagesc(t2,f2,S2);
title('Sound ociscen NO2 optimizacionim metodom');
xlabel('Vreme [s]');
ylabel('Frekvencija [Hz]');
xticks(0:0.2:3);
yticks(0:500:24000);
%Zadavanje parametara za NF3 filtar i odredjivanje prelazne zone
fp13=6000;
fa13=6500;
Fp13=fp13/Fs*2;
Fa13=fa13/Fs*2;
Bt13=(fa13-fp13)/Fs*2*pi;
Bt3=Bt13;
Ap=0.5;
Aa=50;
%Odredjivanje reda filtra
dp=(10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);	
da=10^(-0.05*Aa);
D=(0.005309*log10(dp)*log10(dp)+0.07114*log10(dp)-0.4761)*log10(da);
D=D-(0.00266*log10(dp)*log10(dp)+0.5941*log10(dp)+0.4278);
f=11.01217+0.51244*(log10(dp)-log10(da));
M3=2*pi*D/Bt3-f*Bt3/(2*pi)+1;
M3=ceil(M3);
N3=M3-1;
disp('Pocetni red filtra NF3 je:');
disp(N3);
%Projektovanje NF3 filtra optimizacionim metodom
Hd3=[1 1 0 0];
F3=[0 Fp13 Fa13 1];
h3=firpm(N3,F3,Hd3);
%Crtanje amplitudske karakteristike filtra koji mozda ne zadovoljava zadate
%gabarite
[H3,w3]=freqz(h3,1,256*32);
Ha3=abs(H3);
figure
plot(w3,Ha3),title('Amplitudska karakteristika filtra NF3(firpm)'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
hold on
%Crtanje gabarita
x1h=[0 Fp13*pi]; y1h=[1-dp 1-dp];
x2h=[0 Fp13*pi]; y2h=[1+dp 1+dp];
x1v=[Fp13*pi Fp13*pi]; y1v=[1+dp 0];
x3h=[Fa13*pi pi]; y3h=[da da];
x3v=[Fa13*pi Fa13*pi]; y3v=[da 1];
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r');
kp3=ceil(256*32*Fp13)+1;
ka3=floor(256*32*Fa13)+1;
Hn=Ha3(ka3:256*32)';
Hp=Ha3(1:kp3)';
%Provera da li su zadovoljeni gabariti i ukoliko nisu povecavanje reda
%filtra sve dok ne budu
if((max(Hn)<=da)&&(min(Hp)>=(1-dp))&&(max(Hp)<=(1+dp)))

else
    while(1)
        N3=N3+1;
        M3=M3+1;
        h3=firpm(N3,F3,Hd3);
        [H3,w3]=freqz(h3,1,256*32); 
        Ha3=abs(H3);
        Hn=Ha3(ka3:256*32)';
        Hp=Ha3(1:kp3)';
        if((max(Hn)<=da)&&(min(Hp)>=(1-dp))&&(max(Hp)<=(1+dp)))
            break;
        end
    end
end
disp('Konacni red filtra NF3 je:');
disp(N3);
%Crtanje konacne amplitudkse karakteristike NF3 filtra dobijenog
%optimizacionom metodom
figure
plot(w3,Ha3),title('Amplitudska karakteristika ispravljenog filtra NF3(firpm)'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
hold on
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r');
%Konacno ciscenje zvuka
y3=conv(y2,h3);
%Crtanje spektrograma
[s3,f3,t3]=spectrogram(y3,window,nooverlap,256*32,Fs);
S3=20*log(abs(s3));
figure
imagesc(t3,f3,S3);
title('Sound ociscen optimizacionim metodom');
xlabel('Vreme [s]');
ylabel('Frekvencija [Hz]');
xticks(0:0.2:3);
yticks(0:500:24000);
%Pojacavanje zvuka da bi se lakse cuo i njegovo zapisivanje u novi fajl
y3=2*y3;
audiowrite('sound.wav',y3,Fs);
sound(y3,Fs);
pause
%Zadavanje parametara za NO1 filtar metodom ogranicavanja impulsnog odziva
%i odredjivanje prelazne zone
fp11=1100;
fa11=1300;
fa21=1700;
fp21=1900;
Fp11=fp11/Fs*2;
Fa11=fa11/Fs*2;
Fa21=fa21/Fs*2;
Fp21=fp21/Fs*2;
Bt11=(fa11-fp11)/Fs*2*pi;
Bt21=(fp21-fa21)/Fs*2*pi;
Bt1=min(Bt11,Bt21);
Wc1=fp11/Fs*2*pi+Bt1/2;
Wc2=fp21/Fs*2*pi-Bt1/2;
Ap=0.5;
Aa=40;
Aav=Aa;
Apv=Ap;
%Odredjivanje reda filtra, pravljenje kajzerovog prozora i provera gabarita
uslov1=0;
uslov0=0;
dp1=(10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);
da1=10^(-0.05*Aa);
while(1)
    dp=(10^(0.05*Apv)-1)/(10^(0.05*Apv)+1);
    da=10^(-0.05*Aav);
    d=min(dp,da);
    if d~=da 
        Aav=-20*log10(d);
    end
    beta=0;
    if(Aav>=21 && Aav<=50) 
        beta=0.5842*(Aav-21)^0.4+0.07886*(Aav-21);
    end
    if Aav >50 
        beta=0.1102*(Aav-8.7);
    end
    D=0.9222;
    if Aav>21
        D=(Aav-7.95)/14.36; 
    end
    M=ceil(2*pi*D/Bt1+1);
    N=M-1;  
    prozor=kaiser(M,beta)';

    n=-(M-1)/2:(M-1)/2;
    b=(sin(n*pi)+sin(n*Wc1)-sin(n*Wc2))./(n*pi);
    if (mod(M,2)==1)
        indeks=(M+1)/2 ;
        b(indeks)=1+(Wc1-Wc2)/pi;
    end
    h1=b.*prozor;
    [H1,w1]=freqz(h1,1,256*32);
    Ha1=abs(H1);
    kp11=ceil(256*32*Fp11)+1;
    ka11=floor(256*32*Fa11)+1;
    ka21=ceil(256*32*Fa21)+1;
    kp21=floor(256*32*Fp21)+1;
    Hn=Ha1(ka11:ka21)';
    Hp=[Ha1(1:kp11)' Ha1(kp21:256*32)'];
    %Provera gabarita
    if(max(Hn)<=da1)
        uslov1=1;
    else
        Aav=Aav+0.1;
    end
    if ((min(Hp)>=(1-dp1))&&(min(Hp)>=(1-dp1))&&(max(Hp)<=(1+dp1)))
        uslov2=1;
    else
        Apv=Apv-0.1;
    end   
    if((uslov1==1)&&(uslov2==1))
        disp('Konacni red filtra NO1 je:');
        disp(N);
        break
    end
end
%Cranje amplitude konacnog NO1 filtra imp metodom sa gabaritima
figure
plot(w1,Ha1),title('Amplitudska karakteristika filtra NO1'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
hold on
x1h=[0 Fp11*pi]; y1h=[1-dp1 1-dp1]; 
x2h=[0 Fp11*pi]; y2h=[1+dp1 1+dp1];
x1v=[Fp11*pi Fp11*pi]; y1v=[1+dp1 0];
x3h=[Fa11*pi Fa21*pi]; y3h=[da1 da1];
x3v=[Fa11*pi Fa11*pi]; y3v=[da1 1+dp1];
x4v=[Fa21*pi Fa21*pi]; y4v=[da1 1+dp1];
x4h=[Fp21*pi pi]; y4h=[1+dp1 1+dp1];
x5h=[Fp21*pi pi]; y5h=[1-dp1 1-dp1];
x5v=[Fp21*pi Fp21*pi]; y5v=[0 1+dp1];
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r',x4v,y4v,'r',x4h,y4h,'r',x5h,y5h,'r',x5v,y5v,'r');
figure
stem(0:M-1,h1);
title('Impulsni odziv filtra NO1 dobijen imp metodom');
%Prvo ciscenje zvuka filtrom NO1 dobijenim imp metodom
y1=conv(corrupt',h1);
%Crtanje spektrograma novodobijenog zvuka
[s1,f1,t1]=spectrogram(y1,window,nooverlap,256*32,Fs);
S1=20*log(abs(s1));
figure
imagesc(t1,f1,S1);
title('Sound ociscen NO1 filtrom imp metodom');
xlabel('Vreme [s]');
ylabel('Frekvencija [Hz]');
xticks(0:0.2:3);
yticks(0:500:24000);
%Zadavanje parametara za dobijanje NO2 filtra imp metodom i odredjivanje
%prelazne zone
fp12=2200;
fa12=2400;
fa22=2900;
fp22=3100;
Fp12=fp12/Fs*2;
Fa12=fa12/Fs*2;
Fa22=fa22/Fs*2;
Fp22=fp22/Fs*2;
Bt12=(fa12-fp12)/Fs*2*pi;
Bt22=(fp22-fa22)/Fs*2*pi;
Bt2=min(Bt12,Bt22);
Wc1=fp12/Fs*2*pi+Bt2/2;
Wc2=fp22/Fs*2*pi-Bt2/2;
Ap=0.5;
Aa=40;
%Odredjivanje reda filtra, kajzerovog prozora i provera da li su gabariti
%zadovoljeni
Aav=Aa;
Apv=Ap;
uslov1=0;
uslov0=0;
dp2=(10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);
da2=10^(-0.05*Aa);
while(1)
    dp=(10^(0.05*Apv)-1)/(10^(0.05*Apv)+1);
    da=10^(-0.05*Aav);
    d=min(dp,da);
    if d~=da 
        Aav=-20*log10(d);
    end
    beta=0;
    if(Aav>=21 && Aav<=50) 
        beta=0.5842*(Aav-21)^0.4+0.07886*(Aav-21);
    end
    if Aav >50 
        beta=0.1102*(Aav-8.7);
    end
    D=0.9222;
    if Aav>21
        D=(Aav-7.95)/14.36; 
    end
    M=ceil(2*pi*D/Bt2+1);
    N=M-1;  
    prozor=kaiser(M,beta)';
    n=-(M-1)/2:(M-1)/2;
    b=(sin(n*pi)+sin(n*Wc1)-sin(n*Wc2))./(n*pi);
    if (mod(M,2)==1)
        indeks=(M+1)/2 ;
        b(indeks)=1+(Wc1-Wc2)/pi;
    end
    h2=b.*prozor;
    [H2,w2]=freqz(h2,1,256*32);
    Ha2=abs(H2);
    kp12=ceil(256*32*Fp12)+1;
    ka12=floor(256*32*Fa12)+1;
    ka22=ceil(256*32*Fa22)+1;
    kp22=floor(256*32*Fp22)+1;
    Hn=Ha2(ka12:ka22)';
    Hp=[Ha2(1:kp12)' Ha1(kp22:256*32)'];
    %Provera gabarita
    if(max(Hn)<=da2)
        uslov1=1;
    else
        Aav=Aav+0.1;
    end
    if ((min(Hp)>=(1-dp2))&&(min(Hp)>=(1-dp2))&&(max(Hp)<=(1+dp2)))
        uslov2=1;
    else
        Apv=Apv-0.1;
    end   
    if((uslov1==1)&&(uslov2==1))
        disp('Konacni red filtra NO2 je:');
        disp(N);
        break
    end
end
%Crtanje amplitudkse karakteristike konacnog filtra NO2 dobijenim imp
%metodom sa gabaritima
figure
plot(w2,Ha2),title('Amplitudska karakteristika filtra NO2'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
hold on
x1h=[0 Fp12*pi]; y1h=[1-dp2 1-dp2]; 
x2h=[0 Fp12*pi]; y2h=[1+dp2 1+dp2];
x1v=[Fp12*pi Fp12*pi]; y1v=[1+dp2 0];
x3h=[Fa12*pi Fa22*pi]; y3h=[da2 da2];
x3v=[Fa12*pi Fa12*pi]; y3v=[da2 1+dp2];
x4v=[Fa22*pi Fa22*pi]; y4v=[da2 1+dp2];
x4h=[Fp22*pi pi]; y4h=[1+dp2 1+dp2];
x5h=[Fp22*pi pi]; y5h=[1-dp2 1-dp2];
x5v=[Fp22*pi Fp22*pi]; y5v=[0 1+dp2];
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r',x4v,y4v,'r',x4h,y4h,'r',x5h,y5h,'r',x5v,y5v,'r');
figure
stem(0:M-1,h2);
title('Impulsni odziv filtra NO2 dobijen imp metodom');
%Drugo ciscenje signala filtrom NO2 dobijenim imp metodom
y2=conv(y1,h2);
%Crtanje spektrograma novodobijenog zvuka
[s2,f2,t2]=spectrogram(y2,window,nooverlap,256*32,Fs);
S2=20*log(abs(s2));
figure
imagesc(t2,f2,S2);
title('Sound ociscen NO2 filtrom imp metodom');
xlabel('Vreme [s]');
ylabel('Frekvencija [Hz]');
xticks(0:0.2:3);
yticks(0:500:24000);
%Zadavanje parametara i prelazne zone NF3 filtra dobijenim imp metodom i 
fp13=6000;
fa13=6500;
Fp13=fp13/Fs*2;
Fa13=fa13/Fs*2;
Bt13=(fa13-fp13)/Fs*2*pi;
Bt3=Bt13;
Wc=(2*pi*fp13/Fs+2*pi*fa13/Fs)/2;
Ap=0.5;
Aa=50;
%Odredjivanje reda filtra, pravljenje kajzerovog prozora i provera gabarita
Aav=Aa;
Apv=Ap;
uslov1=0;
uslov0=0;
dp3=(10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);
da3=10^(-0.05*Aa);
while(1)
    dp=(10^(0.05*Apv)-1)/(10^(0.05*Apv)+1);
    da=10^(-0.05*Aav);
    d=min(dp,da);
    if d~=da 
        Aav=-20*log10(d);
    end
    beta=0;
    if(Aav>=21 && Aav<=50) 
        beta=0.5842*(Aav-21)^0.4+0.07886*(Aav-21);
    end
    if Aav >50 
        beta=0.1102*(Aav-8.7);
    end
    D=0.9222;
    if Aav>21
        D=(Aav-7.95)/14.36; 
    end
    M=ceil(2*pi*D/Bt3+1);
    N=M-1;  
    prozor=kaiser(M,beta)';
    n=-(M-1)/2:(M-1)/2;
    b=(sin(n*Wc))./(n*pi);
    if (mod(M,2)==1)
        indeks=(M+1)/2 ;
        b(indeks)=Wc/pi;
    end
    h3=b.*prozor;
    [H3,w3]=freqz(h3,1,256*32);
    Ha3=abs(H3);
    kp3=ceil(256*32*Fp13)+1;
    ka3=floor(256*32*Fa13)+1;
    Hn=Ha3(ka3:256*32)';
    Hp=Ha3(1:kp3)';
    %Provera gabarita
    if(max(Hn)<=da3)
        uslov1=1;
    else
        Aav=Aav+0.1;
    end
    if ((min(Hp)>=(1-dp3))&&(min(Hp)>=(1-dp3))&&(max(Hp)<=(1+dp3)))
        uslov2=1;
    else
        Apv=Apv-0.1;
    end   
    if((uslov1==1)&&(uslov2==1))
        disp('Konacni red filtra NF3 je:');
        disp(N);
        break
    end
end
%Crtanje amplitude konacnog filtra NF3 dobijenog imp metodom sa gabaritima
figure
plot(w3,Ha3),title('Amplitudska karakteristika filtra NF3'),grid on;
xlabel('w');
ylabel('|H(e(jw))|');
hold on
x1h=[0 Fp13*pi]; y1h=[1-dp3 1-dp3];
x2h=[0 Fp13*pi]; y2h=[1+dp3 1+dp3];
x1v=[Fp13*pi Fp13*pi]; y1v=[1+dp3 0];
x3h=[Fa13*pi pi]; y3h=[da3 da3];
x3v=[Fa13*pi Fa13*pi]; y3v=[da3 1];
plot(x1h,y1h,'r',x1v,y1v,'r',x2h,y2h,'r',x3h,y3h,'r',x3v,y3v,'r');
figure
stem(0:M-1,h3);
title('Impulsni odziv filtra NF3 dobijen imp metodom');
%Konacno ciscenje zvuka
y3=conv(y2,h3);
%Crtanje spektrograma konacnog zvuka
[s3,f3,t3]=spectrogram(y3,window,nooverlap,256*32,Fs);
S3=20*log(abs(s3));
figure
imagesc(t3,f3,S3);
title('Sound ociscen filtrom NO3 imp metodom (konacan)');
xlabel('Vreme [s]');
ylabel('Frekvencija [Hz]');
xticks(0:0.2:3);
yticks(0:500:24000);
%Pojacavanje malo zvuka da bi se lepse cuo i zapisivanje u novi fajl
y4=2*y3;
audiowrite('out_signal_2017_316.wav',y4,Fs);
sound(y4,Fs);