function newMAT=isi_shuffle(MAT)
newMAT=zeros(size(MAT));
for i=1:size(MAT,1)
    ix=  find(MAT(i,:));
    if ~isempty(ix)
        dix= [ix(1) diff(ix)];
        [~,rix]= sort(rand(1,length(dix)));
        ix=dix(rix);
        new_ix= cumsum([ix]);
        newMAT(i,new_ix)=1;
    end
end
end