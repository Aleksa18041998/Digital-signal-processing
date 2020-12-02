function h=power_line_noise_filter(fs,fc,Aa,Ap)
    %Odredjivanje parametara analognog prototipa i predistorzija 
    T=1/fs;
    Apv=Ap;
    Aav=Aa;
    fp1=fc-2;
    fp2=fc+2;
    fa1=fc-0.5;
    fa2=fc+0.5;
    Wp1=2*pi*fp1;
    Wp2=2*pi*fp2;
    Wa1=2*pi*fa1;
    Wa2=2*pi*fa2;
    B=Wp2-Wp1;
    W0=sqrt(Wp1*Wp2);
    WaN=B*Wa1/(W0^2-Wa1^2);
    Wp1=2/T*tan(Wp1/2/fs);
    Wp2=2/T*tan(Wp2/2/fs);
    Wa1=2/T*tan(Wp1/2/fs);
    Wa2=2/T*tan(Wa2/2/fs);
    B=2/T*tan(B/2/fs);
    W0=2/T*tan(W0/2/fs);  
    WaN=2/T*tan(WaN/2/fs);
    WpN=1;
    %Odredjivanje reda filtra i dobijanje analognog prototipa, digitalnog
    %filra
    while(1)
        k=sqrt(1-(WpN/WaN)^2);
        D=(10^(0.1*Aav)-1)/(10^(0.1*Apv)-1);
        q0=(1/2)*((1-sqrt(k))/(1+sqrt(k)));
        q=q0+2*q0^5+15*q0^9+15*q0^13;
        N=ceil(log10(16*D)/(log10(1/q)));
        [zn,pn,kn]=ellipap(N,Apv,Aav);
        baN=kn*poly(zn);
        aaN=poly(pn);
        [b,a]=lp2bs(baN,aaN,W0,B);
        [bd,ad]=bilinear(b,a,fs);
        h=[bd' ad'];
        Nfreqz=20000;
        [hd,Wd]=freqz(bd,ad,Nfreqz);
        Hd=abs(hd);
        [phi,Wd1]=phasez(bd,ad,Nfreqz);
        fd=(Wd*fs)/(2*pi);
        %provera da li su uspesno zadovoljeni gabariti
        uslov1=0;
        uslov2=0;
        df=(fs/2)/Nfreqz;
        ip1=floor(fp1/df)+1;
        ia1=ceil(fa1/df)+1;
        ia2=ceil(fa2/df)+1;
        ip2=floor(fp2/df)+1;
        Ha=Hd(ia1:ia2);
        Hp=[Hd(1:ip1)' Hd(ip2:length(Hd))'];
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
    %Crtanje frekvencijskih karakteristika
    f=Wd/(2*pi)*fs;
    figure
    plot(f,20*log10(Hd),'LineWidth',2);
    title('Amplitudska karakteristika NO filtra u logaritamskoj razmeri'),grid on;
    xlabel('Ucestanost (Hz)');
    ylabel('|H(z)|');
    %Crtanje zadatih gabarita 
    hold on
    x1h=[0 fp1]; y1h=[-Ap -Ap];
    x1v=[fp1 fp1]; y1v=[0 -Ap];
    x2h=[fa1 fa2]; y2h=[-Aa -Aa];
    x21v=[fa1 fa1]; y21v=[-Aa min(20*log10(Hd))];
    x22v=[fa2 fa2]; y22v=[-Aa min(20*log10(Hd))];
    x3h=[fp2 fs/2]; y3h=[-Ap -Ap];
    x3v=[fp2 fp2]; y3v=[0 -Ap];
    plot(x1h,y1h,'r',x2h,y2h,'r',x1v,y1v,'r',x21v,y21v,'r',x22v,y22v,'r',x3h,y3h,'r',x3v,y3v,'r','LineWidth',1.5);
    hold off
    %Crtanje nula i polova filtra
    figure
    zo=roots(bd);
    zp=roots(ad);
    [hz,hp,ht]=zplane(zo,zp); 
    set(findobj(hp,'Type','line'),'LineWidth',3,'MarkerSize',12);
    title('Raspored polova projektovanog NO filtra')
    xlabel('Re(s)');
    ylabel('Im(s)');
    f1=Wd1/(2*pi)*fs;
    %Crtanje fazne karakteristike
    figure
    plot(f1,phi),title('Fazna karakteristika NO filtra'),grid on;

    
    
    