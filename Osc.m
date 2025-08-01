function AD=Osc(imf,kbest,fs,xuhao)
for i=1:kbest
    [p,tp]=regularity_new(imf(:,i),fs,0);
    AD(i,:)=[xuhao(i) tp p];
end
end