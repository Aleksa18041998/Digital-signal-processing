function    y=FI_IIR_direct_II_cascade(b,a,x)
    
    fi_params=struct('pvOUT_SIGNAL_BITLENGTH',x.WordLength+a.WordLength,'pvOUT_SIGNAL_FRAC',x.FractionLength+a.FractionLength,...
    'yOUT_SIGNAL_BITLENGTH',x.WordLength+a.WordLength+b.WordLength,'yOUT_SIGNAL_FRAC',x.FractionLength+a.FractionLength+b.FractionLength);
    FixedPointAttributes=fimath('RoundingMethod','Floor','OverflowAction','Wrap',...
    'ProductMode','SpecifyPrecision','ProductWordLength',fi_params.yOUT_SIGNAL_BITLENGTH ,'ProductFractionLength',fi_params.yOUT_SIGNAL_FRAC,...
    'SumMode','SpecifyPrecision','SumWordLength',fi_params.yOUT_SIGNAL_BITLENGTH,'SumFractionLength',fi_params.yOUT_SIGNAL_FRAC);

    y=fi(zeros(1,length(x)),true,fi_params.yOUT_SIGNAL_BITLENGTH,fi_params.yOUT_SIGNAL_FRAC,FixedPointAttributes);
    p=fi(zeros(1,length(x)),true,fi_params.pvOUT_SIGNAL_BITLENGTH,fi_params.pvOUT_SIGNAL_FRAC,FixedPointAttributes);  
    v=fi(zeros(1,length(x)),true,fi_params.pvOUT_SIGNAL_BITLENGTH,fi_params.pvOUT_SIGNAL_FRAC,FixedPointAttributes);
    delay_line=fi(zeros(length(a),1),x.Signed,x.WordLength,x.FractionLength,FixedPointAttributes); 
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