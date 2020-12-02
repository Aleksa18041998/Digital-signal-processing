function y=IIR_direct_II_cascade(b,a,x)
    [brvrsta,brkolona]=size(b);
    ncascade=brvrsta;
    for i=1:ncascade
        delay_line=zeros(brkolona-1,1);
        v(1)=x(1);
        for k=2:length(x)
            delay_line=[v(k-1);delay_line(1:end-1)];
            p(k)=-a(i,:)*delay_line;
            v(k)=p(k)+x(k);
        end
        delay_line=zeros(brkolona,1);
        for k=1:length(x)
            delay_line=[v(k);delay_line(1:end-1)];
            y(k)=b(i,:)*delay_line;
            x(k)=y(k);
        end
    end