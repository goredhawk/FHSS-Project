function [ IM2,LINE ] = improcess( IM, SE)
%形态学滤波处理
%   IM：输入的二值图像
%   SE：形态学滤波的结构元素
%   IM2：滤波后的图像
%   LINE：检测到的直线坐标

IM2=IM;
for i=1:2
    %开运算
    IM2=imdilate(IM2,SE);
    IM2=imerode(IM2,SE);
    %闭运算
    IM2=imerode(IM2,SE);
    IM2=imdilate(IM2,SE);
end

[h,w]=size(IM2);
LINE=zeros(1,w);
for j=1:w
    k=0;
    LINE(j)=0;
    for i=1:h
        if IM2(i,j) == 1
            LINE(j)=LINE(j)+i;
            k=k+1;
        end
    end
    if k ~=0
        LINE(j)=LINE(j)/k;
    end
end

end

