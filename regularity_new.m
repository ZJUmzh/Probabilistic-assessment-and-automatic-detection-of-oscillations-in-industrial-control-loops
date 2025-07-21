function [p,aTp]=regularity_new(x,fs,flag)
if std(x)>0
n=length(x);
acf_num=acf(x,n-1);
zeroCrossings = diff(sign(acf_num)) ~= 0;
zeroCrossings = find(zeroCrossings) + 1;
d=1:n-1;
q=(1.96)*(1./sqrt(n-d));
nmax=0;
for i=1:n-1
    if abs(acf_num(i))>q(i) 
        nmax=i;
    end
end
if max(zeroCrossings)<nmax
    zeroCrossings=zeroCrossings;
else
    nmin=min(find(zeroCrossings>nmax));
    zeroCrossings=zeroCrossings((1:nmin-1),1);
end
if length(zeroCrossings)>=3 
    Tp=2*diff(zeroCrossings)/fs;
    aTp=mean(Tp);
    std_Tp=std(Tp);
    L=length(Tp);
 p=probability(L,std_Tp,aTp);
else
p=NaN;
aTp=NaN;
end
if flag==1
figure;
plot(acf_num,'linewidth',2,'color', '#2679F4');
hold on;
plot(q,'Color', 'r');
hold on;
plot(-q,'Color', 'r');
hold on;
if isempty(zeroCrossings)==0
plot(zeroCrossings,0,'k.','MarkerSize',15);
end
set(gca,'FontName','微软雅黑','FontSize',24);
set(gca,'box','on','linewidth',1.5);
set(gcf,'unit','centimeters','position',[0 0 25 18]);
xlabel('lags','FontSize',24,'FontName','微软雅黑','color', '#333333');
ylabel('ACF','FontSize',24,'FontName','微软雅黑','color', '#333333');
end
else
aTp=NaN;
p=NaN;
end
end