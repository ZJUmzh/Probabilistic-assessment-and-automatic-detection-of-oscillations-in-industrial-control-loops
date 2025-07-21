function [AD_TABLE,AD]=Osc(imf,kbest,fs,xuhao)
for i=1:kbest
    [p,tp]=regularity_new(imf(:,i),fs,0);
    AD(i,:)=[xuhao(i) tp p];
end


AD_NAME = {'Modal number';'Oscillation period/s';'Oscillation possibility'};
AD=AD';
AD_TABLE=table(AD_NAME,AD);

end