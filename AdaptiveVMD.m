function [imf,kbest,kmax]=AdaptiveVMD(x,fs)
k=2;
rng(0);
[~,~,info]=vmd(x,'NumIMF' ,k);
max_fs_c=max(info.CentralFrequencies*fs);
if max_fs_c>=fs/2.56
    kmax=1;
else
    while max_fs_c<(fs/2.56) 
        k=k+1;
        [~,~,info]=vmd(x,'NumIMF' ,k);
        max_fs_c=max(info.CentralFrequencies*fs);
    end
    kmax=k-1;
end
if kmax==1 
    imf=x;
    kbest=1;
else
    Es=sum(x.^2); 
    for i=2: kmax 
        [imf,~,info]=vmd(x,'NumIMF' ,i);
        MI_sum=0;
        E=0;
        for j=1:i 
            MI_sum=MI_sum+mutualinformation(imf(:,j),x);
            E=E+sum(imf(:,j).^2);
            MI1(j)=getSampEn(imf(:,j));
        end
        MI_ave2(i)=MI_sum/(i);
        e(i)=abs(Es-E)/Es;
        MI_ave1(i)=min(MI1);
        cv(i)=std(info.CentralFrequencies*fs)/mean(info.CentralFrequencies*fs);
        MI_c(i)=cv(i)*MI_ave2(i)/MI_ave1(i)/e(i);
    end
    kbest=min(find(MI_c==max(MI_c)));
    imf=vmd(x,'NumIMF' ,kbest);
end
end