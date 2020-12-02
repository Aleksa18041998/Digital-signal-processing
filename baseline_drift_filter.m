function h=baseline_drift_filter(fs,fa,fp,Aa,Ap)
    %Odredjivanje parametara filtra i predistorzija
    Wp=2*pi*fp;
    Wa=2*pi*fa;
    T=1/fs;
    Apv=Ap;
    Aav=Aa;
    Wp=2/T*tan(Wp/2/fs);
    Wa=2/T*tan(Wa/2/fs);
    WpN=1;
    WaN=Wp/Wa;
    %Petlja kojom se odredjuje potreban red filtra i provera da li
    %zadovoljava zadate pocetne gabarite
    while(1)
        D=(10^(0.1*Aav)-1)/(10^(0.1*Apv)-1);
        k=WpN/WaN;
        N=ceil(acosh(sqrt(D))/acosh(1/k));
        [zn,pn,kn]=cheb2ap(N,Aav);
        bn=kn*poly(zn);
        an=poly(pn);
        [b,a]=lp2hp(bn,an,Wa);
        [bd,ad]=bilinear(b,a,fs);
        h=[bd' ad'];
        Nfreqz=20000;
        [hd,Wd]=freqz(bd,ad,Nfreqz);
        Hd=abs(hd);
        [phi,Wd1]=phasez(bd,ad,Nfreqz);
        fd=(Wd*fs)/(2*pi);
        %provera da li su zadovoljeni gabariti
        uslov1=0;
        uslov2=0;
        df=(fs/2)/Nfreqz;
        ia=ceil(fa/df)+1;
        ip=floor(fp/df)+1;
        Ha=Hd(1:ia);
        Hp=Hd(ip:length(Hd));
        if(max(20*log10(Ha))>-Aa)
            Aav=Aav+0.1;
        else
            uslov1=1;
        end
        if(min(20*log10(Hp))<-Ap)
            Apv=Apv-0.1;
        else
            uslov2=1;
        end   
        if((uslov1==1)&&(uslov2==1))
            Apv
            Aav
            N
            break
        end
    end
    %Odredjivanje frekvencijskih karakteristika filtra
    f=Wd/(2*pi)*fs;
    figure
    plot(f,20*log10(Hd),'LineWidth',2),grid on ;
    title('Amplitudska karakteristika VF filtra u logaritamskoj razmeri');
    xlabel('Ucestanost (Hz)');
    ylabel('|H(z)|');
    %Crtanje gabarita koje treba da zadovolji
    hold on
    x1h=[0 fa]; y1h=[-Aa -Aa];
    x1v=[fa fa]; y1v=[-Aa min(20*log10(Hd))];
    x2h=[fp fs/2]; y2h=[-Ap -Ap];
    x2v=[fp fp]; y2v=[0 -Ap];
    plot(x1h,y1h,'r',x2h,y2h,'r',x1v,y1v,'r',x2v,y2v,'LineWidth',1.5);
    hold off
    %Odredjivanje polova filtra
    figure
    zo=roots(bd);
    zp=roots(ad);
    [hz,hp,ht]=zplane(zo,zp); 
    set(findobj(hp,'Type','line'),'LineWidth',3,'MarkerSize',12);
    title('Raspored polova projektovanog VF filtra')
    xlabel('Re(s)');
    ylabel('Im(s)');
    figure
    f1=Wd1/(2*pi)*fs;
    plot(f1,phi),title('Fazna karakteristika VF filtra'),grid on;

    