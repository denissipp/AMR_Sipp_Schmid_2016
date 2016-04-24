function y = rd(x,n)
    y=sign(x)*round(abs(x)*10^(-ceil(log10(abs(x)))+n))/10^(-ceil(log10(abs(x)))+n);
end
