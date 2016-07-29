clc;
clear;
warning('off');
%%%%%   �����ź�����ģ�鼰��ʱ����Ҷ�任ģ��  %%%%%%%
%   ����ź�˵�� 
%  signal_out -- ��ʱ����Ҷ�任�����ź�
%   ����ģʽ˵��
% modulation_mode��'none'      --�޵���
%                  '2fsk'      --2FSK����
%                  '4fsk'      --2FSK����
%                  'qpsk'      --QPSK����
%                  'bpsk'      --BPSK����
%                  '16qam'     --16QAM���� 
fset=100*10^6;                %������Ч����
fs=fset*1.25;                 %���źŲ�����(Hz)����Ϊ��Ч�����1.25��
carrier_freq_separation=6.25*10^6;   %�ز�Ƶ�ʼ��(Hz)
hopping_rate=10000;          %��Ƶ����(hop/s)
modulation_mode='16qam';     %����ģʽ
switch_time=10^-6;          %�л�ʱ��
symbol_rate=1*10^6;         %��������(��������)
signal_length=0.001;        %�ź�ʱ�䳤��(s)
SNR=3;                     %�����

%�����ź�����
Nhop=fix(hopping_rate*signal_length);   %%���������
hopping_cycle=1/hopping_rate; %��Ƶ����
symbol_cycle=1/symbol_rate; %��������
sample_per_symbol=fix(fs/symbol_rate);             %ÿһ�����ŵĲ�������
sample_num_per_hop=fix(fs/hopping_rate);           %ÿһ���Ĳ�������
symbol_num_per_hop=fix(symbol_rate/hopping_rate);  %ÿһ���ķ�����
sample_num_switch_time=fix(fs*switch_time);        %�л�ʱ���ڵĲ�������

f=zeros(Nhop,1);
for i=1:Nhop
    f(i)=randi(15,1,1);
end
x=[];
if strcmp(modulation_mode,'none')
    for i=1:Nhop-1
        carrier_num=f(i);
        sim('none_v2');
        x=[x, reshape(signal_out,1,[])];
        t=0:1/fs:switch_time;
        y1=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',0);
        y2=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',90);
        y=y1-y2*1i;
        y=awgn(y,SNR);
        x =[x, y];
    end
    carrier_num=f(Nhop);
    sim('none_v2');
    x=[x, reshape(signal_out,1,[])];
elseif strcmp(modulation_mode,'2fsk')
    freq_separation=2*10^6; %2fsk����Ƶ�ʼ��
    for i=1:Nhop-1
        carrier_num=f(i);
        sim('fsk2_v2');
        x=[x, reshape(signal_out,1,[])];
        if f(i)~=f(i+1)
            t=0:1/fs:switch_time;
            y1=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',0);
            y2=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',90);
%             y1=y1*0;
%             y2=y2*0;
            y=y1-y2*1i;
            y=awgn(y,SNR);
            x =[x, y];
        end
    end
    carrier_num=f(Nhop);
    sim('fsk2_v2');
    x=[x, reshape(signal_out,1,[])];
elseif strcmp(modulation_mode,'bpsk')
    for i=1:Nhop-1
        carrier_num=f(i);
        sim('bpsk_v2');
        x=[x, reshape(signal_out,1,[])];
        t=0:1/fs:switch_time;
        y1=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',0);
        y2=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',90);
        y=y1-y2*1i;
        y=awgn(y,SNR);
        x =[x, y];
    end
    carrier_num=f(Nhop);
    sim('bpsk_v2');
    x=[x, reshape(signal_out,1,[])];
elseif strcmp(modulation_mode,'16qam')
    for i=1:Nhop-1
        carrier_num=f(i);
        sim('qam16');
        x=[x, reshape(signal_out,1,[])];
        t=0:1/fs:switch_time;
        y1=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',0);
        y2=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',90);
        y=y1-y2*1i;
        y=awgn(y,SNR);
        x =[x, y];
    end
    carrier_num=f(Nhop);
    sim('qam16');
    x=[x, reshape(signal_out,1,[])];
 elseif strcmp(modulation_mode,'qpsk')
    for i=1:Nhop-1
        carrier_num=f(i);
        sim('qpsk');
        x=[x, reshape(signal_out,1,[])];
        t=0:1/fs:switch_time;
        y1=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',0);
        y2=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',90);
        y=y1-y2*1i;
        y=awgn(y,SNR);
        x =[x, y];
    end
    carrier_num=f(Nhop);
    sim('qpsk');
    x=[x, reshape(signal_out,1,[])];
end

%��ʱ����Ҷ�任���ò���(������ΪHamming��)
analysis_window_length=30;
overlap=0;
FFT_length=1024;
%%%Hamming��%%%%%%%
hwin=zeros(1,analysis_window_length);
win_hamming=hamming(analysis_window_length);
for i=1:analysis_window_length
    hwin(i)=win_hamming(i);
end

z=zeros(fix(length(x)/analysis_window_length),FFT_length);
pow=zeros(fix(length(x)/analysis_window_length),1);
fe=zeros(fix(length(x)/analysis_window_length),1);
for i=1:fix(length(x)/analysis_window_length)
    hwin_x = x((i-1)*analysis_window_length+1:i*analysis_window_length).*hwin;
    z(i,:)=abs(fftshift(fft(hwin_x,FFT_length)));
    pow(i)=sum(abs(x((i-1)*analysis_window_length+1:i*analysis_window_length)));
    [a,b]=max(z(i,:));
    fe(i)=b;
end

figure(3);mesh(z);view(2);title('�ٲ�ͼ');
z=z';
z=medfilt2(z);                            %��ֵ�˲�
z=filter2(fspecial('disk',3),z);
[a,z_max]=max(z);
figure(1);plot(fe);title('ʱƵͼ');
hold on;plot(z_max,'r');hold off;
freq_threshold=10;
if strcmp(modulation_mode,'2fsk')
    freq_threshold=fix(freq_separation/fs*1024)+10;
elseif strcmp(modulation_mode,'16qam')
    freq_threshold=15;
end
[ hopPoint ] = timeEstimate( z_max, freq_threshold );   %��������
figure(2);plot(pow);title('ʱ��-����');

%ʱ���������
dwell_time=zeros(fix(length(hopPoint)/2),1);
sw_time=zeros(fix((length(hopPoint)-1)/2),1);
k=1;
for i=2:2:length(hopPoint)
    dwell_time(k)=(hopPoint(i)-hopPoint(i-1))*analysis_window_length/fs;
    k=k+1;
end

k=1;
for i=3:2:length(hopPoint)
    sw_time(k)=(hopPoint(i)-hopPoint(i-1))*analysis_window_length/fs;
    k=k+1;
end

%Ƶ�ʲ�������
k=1;
freq_hop=zeros(fix(length(hopPoint)/2),1);       %FFT�����ز�Ƶ��
freq_hop_avg=zeros(fix(length(hopPoint)/2),1);   %ÿһ���ز�Ƶ�ʵĴ��Թ���ֵ
hopPoint2=(hopPoint-0.5)*analysis_window_length;
for i=2:2:length(hopPoint2)
    freq_hop_avg(k)=(mean(z_max(hopPoint(i-1):hopPoint(i)))-FFT_length/2-1)/FFT_length*fs; 
    if strcmp(modulation_mode,'bpsk')
        t=0:1/fs:(hopPoint2(i)-hopPoint2(i-1))/fs;
        y=x(hopPoint2(i-1):hopPoint2(i)).^2;
        y=y.*(cos(4*pi*freq_hop_avg(k)*t)-sin(4*pi*freq_hop_avg(k)*t)*1i);
        fft_hop=abs(fftshift(fft(y,100000)));
        [a,b]=max(fft_hop);
        freq_hop(k)=(b-50001)/100000*fs/2+freq_hop_avg(k);     
    elseif strcmp(modulation_mode,'none')
        fft_hop=abs(fftshift(fft(x(hopPoint2(i-1):hopPoint2(i)),100000)));
        [a,b]=max(fft_hop);
        freq_hop(k)=(b-50001)/100000*fs;
    elseif strcmp(modulation_mode,'16qam')
        t=0:1/fs:(hopPoint2(i)-hopPoint2(i-1))/fs;
        y=x(hopPoint2(i-1):hopPoint2(i));
        y=real(y.*(cos(4*pi*freq_hop_avg(k)*t)-sin(4*pi*freq_hop_avg(k)*t)*1i)).^2;
        fft_hop=abs(fftshift(fft(y,100000)));
        [a,b]=max(fft_hop);
        if freq_hop_avg(k)<0
            freq_hop(k)=(b-50001)/100000*fs;
        else
            freq_hop(k)=-(b-50001)/100000*fs;
        end
    elseif strcmp(modulation_mode,'qpsk')
        t=0:1/fs:(hopPoint2(i)-hopPoint2(i-1))/fs;
        y=x(hopPoint2(i-1):hopPoint2(i)).^4;
        y=y.*(cos(8*pi*freq_hop_avg(k)*t)-sin(8*pi*freq_hop_avg(k)*t)*1i);
        fft_hop=abs(fftshift(fft(y,100000)));
        [a,b]=max(fft_hop);
        freq_hop(k)=(b-50001)/100000*fs/4+freq_hop_avg(k);
    elseif strcmp(modulation_mode,'2fsk')
        fft_hop=abs(fftshift(fft(x(hopPoint2(i-1):hopPoint2(i)),100000)));
        [a,b1]=max(fft_hop);
        b=fix(freq_separation/fs*100000);      %Ƶ�ʼ����Ӧ�ĵ���
        if b1-b < 0
            [a22,b22]=max(fft_hop(b1+0.5*b:b1+1.5*b));
            b2=b1+b22+0.5*b-1;
        elseif b1+b >100000
            [a21,b21]=max(fft_hop(b1-1.5*b:b1-0.5*b));
            b2=b1+b21-1.5*b-1;
        else
            [a21,b21]=max(fft_hop(b1-1.5*b:b1-0.5*b));
            [a22,b22]=max(fft_hop(b1+0.5*b:b1+1.5*b));
            if a21 > a22
                b2=b1+b21-1.5*b-1;
            else
                b2=b1+b22+0.5*b-1;
            end
        end
        freq_hop(k)=((b1+b2)/2-50001)/100000*fs;
    end
    k=k+1;
end

% freq_hop1=zeros(fix(length(hopPoint)/2),1);
% %���ֵ����źŵ��ز����ƣ�kay��
% for i=2:2:length(hopPoint2)
%     angle_time=angle(x(hopPoint2(i-1)+60:hopPoint2(i)-60)); %˲ʱ��λ
%     freq_time=zeros(length(angle_time)-2,1);
%     for j=2:length(angle_time)
%         if abs(angle_time(j)-angle_time(j-1)) < pi
%             freq_time(j-1)=angle_time(j)-angle_time(j-1);
%         else
%             if (angle_time(j)-angle_time(j-1)) > 0
%                 freq_time(j-1)=angle_time(j)-angle_time(j-1)-2*pi;
%             else
%                 freq_time(j-1)=angle_time(j)-angle_time(j-1)+2*pi;
%             end
%         end
%     end
% %     freq_time=medfilt1(freq_time,3);
%     N=length(freq_time);
%     mean_freq=mean(freq_time);
%     sum=0;
%     k=0;
%     for j=1:N
%         if abs(freq_time(j)-mean_freq)< pi/2
%             sum=sum+freq_time(j);
%             k=k+1;
%         end
%     end
% %     freq_time=freq_time.*((1-(2*k-N)).^2/N^2)*3/2/N;
%     freq_hop1=sum/k*fs/2/pi
% end


f=(f-7.5)*6.25*10^6;





%����ģ��
% k=1;
% dehop_y=[];
% for i=2:2:length(hopPoint2)
%     t=0:1/fs:(hopPoint2(i)-hopPoint2(i-1))/fs;
%     dehop_x=x(hopPoint2(i-1):hopPoint2(i)).*(cos(2*pi*freq_hop(k)*t)-sin(2*pi*freq_hop(k)*t)*1i);
%     dehop_y=[dehop_y,dehop_x];
%     k=k+1;
% end





    