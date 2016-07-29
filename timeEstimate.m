function [ hopPoint ] = timeEstimate( FE, THRESHOLD )
%UNTITLED2 Summary of this function goes here
%   FE:载波频率对应点
%   THRESHOLD:跳变门限

%预处理
for i=2:length(FE)-1
    if abs(FE(i)-FE(i-1)) > THRESHOLD
        if abs(FE(i+1)-FE(i-1)) < THRESHOLD
            FE(i)=FE(i-1);
        end
    end        
end

%第一跳
for i=1:length(FE)-19
    k=0;
    for j=1:19
        if FE(i+j)-FE(i) < THRESHOLD
            k=k+1;
        else break;
        end
    end
    if k==19
        hopPoint(1)=i;
        break;
    end
end

%跳变点
k=2;
n=0;
flag=1;%1--寻找跳变点   0--寻找起始点
if i<length(FE)-19
    ref_freq=mean(FE(hopPoint(1):hopPoint(1)+19));
    for j=hopPoint(1)+20:length(FE)-19
        if 1==flag
            if abs(FE(j)-ref_freq)>THRESHOLD && abs(FE(j+1)-ref_freq) > THRESHOLD && abs(FE(j+2)-ref_freq) > THRESHOLD
                hopPoint(k)=j;
                k=k+1;
                flag=0;
            end
        else
            for jj=0:18
                if abs(FE(j+jj)-FE(j-1)) < THRESHOLD
                    n=n+1;
                else break;
                end
                if 19==n
                    hopPoint(k)=j-1;
                    ref_freq=mean(FE(hopPoint(k):hopPoint(k)+19));
                    k=k+1;
                    flag=1;
                end
            end
            n=0;
        end
    end
else
    hopPoint=[];
end
end

