function [ y ] = stft_get_carrier( fft_value )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[max_value,max_value_index]=max(fft_value);
len=length(fft_value);
sum_value=sum(fft_value);
new_index=1-max_value_index:len-max_value_index;
for i=1:len
    if abs(new_index(i)) > len/2
        if new_index(i) >0
            new_index(i)=new_index(i)-len;
        else
            new_index(i)=new_index(i)+len;
        end
    end
end
fft_value=fft_value.*new_index;
y=sum(fft_value)/sum_value+max_value_index;
end

