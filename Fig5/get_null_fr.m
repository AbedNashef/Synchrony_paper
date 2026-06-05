function nullfr= get_null_fr(realfr,bin)
ix= find(realfr);
nullfr=zeros(1,length(realfr));
bin=1;
for i=2:length(ix)-1
    T1=ix(i-1);
    T=ix(i);
    T2=ix(i+1);
    for t=round(T1+.5*(T-T1)):round(T+.5*(T2-T))
        nullfr(t)=(2*bin)/(T2-T1);
    end
end

end