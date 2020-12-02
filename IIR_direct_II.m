function y=IIR_direct_II(b,a,x)
    delay_line=zeros(length(a),1);
    v(1)=x(1);
    for k=2:length(x)
        delay_line=[v(k-1);delay_line(1:end-1)];
        p(k)=-a*delay_line;
        v(k)=p(k)+x(k);
    end
    delay_line=zeros(length(b),1);
    for k=1:length(x)
        delay_line=[v(k);delay_line(1:end-1)];
        y(k)=b*delay_line;
    end