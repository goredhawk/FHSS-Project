function [ y ] = kay( x,fs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% len=length(x);

angle_x=angle(conj(x(1:end-1)).*x(2:end));
angle_x_avg=mean(angle_x);
sum=0;
k=0;
for i=1:length(angle_x)
    if abs(angle_x(i)-angle_x_avg)< pi/2
        sum=sum+angle_x(i);
        k=k+1;
    end
end
y=sum/k/pi*1024;
end

