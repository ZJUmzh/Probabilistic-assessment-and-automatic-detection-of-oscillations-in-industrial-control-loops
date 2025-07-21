function [imf_new,xuhao,n,delta]=significant_IMFs(imf,x,kbest)
if kbest>1
    for i=1:kbest
        rr=corrcoef(imf(:,i),x);
        r(i,1)=rr(1,2);
        r(i,2)=i;
    end
    r=sortrows(r,1,'descend');
    for i=1:kbest
        if i==1
            RR(i)=r(i,1);
            imf_leiji=imf(:,r(i,2));
        else
            imf_leiji=imf_leiji+imf(:,r(i,2));
            TT=corrcoef(imf_leiji,x);
            RR(i)=TT(1,2);
        end
    end
    n=0;
    for i=1:kbest
        if i==1
            delta(1)=1;
        else
            delta(i)=(RR(i)-RR(i-1))/RR(i-1);
        end
        if delta(i)>=0.05
            n=n+1;
            imf_new(:,n)=imf(:,r(i,2));
            xuhao(n)=r(i,2);
        end
    end
    figure
    edges = [1:1:length(delta)];
    yyaxis left;bar(edges,delta,1,'EdgeColor','k','Linewidth',1,'FaceColor',[0.2 0.6 0.8]);
    hold on;
    plot([0,length(delta)+1],[0.05,0.05],'color', [0 0.447 0.741] ,'linewidth',3)
    xlim([0,length(delta)+1]);
    hold on;
    yyaxis right;plot(RR,'-o','linewidth',2);
    set(gca,'FontName','微软雅黑','FontSize',24);
    set(gca,'box','on','linewidth',1.5);
    set(gcf,'unit','centimeters','position',[0 0 25 18]);
else
    imf_new=imf;
    xuhao=1;
    n=1;
    delta=1;
    RR=0;
    r=0;
end
end
