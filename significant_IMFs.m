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
else
    imf_new=imf;
    xuhao=1;
    n=1;
    delta=1;
    RR=0;
    r=0;
end
end
