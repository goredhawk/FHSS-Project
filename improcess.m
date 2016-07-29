function [ IM2,LINE ] = improcess( IM, SE)
%��̬ѧ�˲�����
%   IM������Ķ�ֵͼ��
%   SE����̬ѧ�˲��ĽṹԪ��
%   IM2���˲����ͼ��
%   LINE����⵽��ֱ������

IM2=IM;
for i=1:2
    %������
    IM2=imdilate(IM2,SE);
    IM2=imerode(IM2,SE);
    %������
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

