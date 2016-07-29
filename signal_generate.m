function y=signal_generate(hopping_rate, switch_time, carrier_freq_separation, modulation_mode, symbol_rate, SNR, signal_length, fs)
% Frequency-Hopping Signals Generater
% hopping_rate:                  ��Ƶ����(hop/s)
% switch_time:                   �л�ʱ��
% carrier_freq_separation:       �ز�Ƶ�ʼ��(Hz)
% modulation_mode:               ����ģʽ
% symbol_rate:                   ��������
% SNR:                           �����
% signal_length:                 �ź�ʱ�䳤��(s)
% fs:                            ���źŲ�����(Hz)

Nhop=fix(hopping_rate*signal_length);   %%���������
hopping_cycle=1/hopping_rate; %��Ƶ����
symbol_cycle=1/symbol_rate; %��������
sample_per_symbol=fix(fs/symbol_rate);             %ÿһ�����ŵĲ�������
sample_num_per_hop=fix(fs/hopping_rate);           %ÿһ���Ĳ�������
symbol_num_per_hop=fix(symbol_rate/hopping_rate);  %ÿһ���ķ�����
sample_num_switch_time=fix(fs*switch_time);

f=zeros(Nhop,1);
for i=1:Nhop
    f(i)=randi(15,1,1);
end
y=[];
if strcmp(modulation_mode,'none')
    
elseif strcmp(modulation_mode,'2fsk')
    freq_separation=0.5*10^6; %2fsk����Ƶ�ʼ��
    for i=1:Nhop-1
        carrier_num=f(i);
        sim('fsk2_v2');
        y=[y, reshape(signal_out,1,[])];
        t=0:1/fs:switch_time;
        y1=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',0);
        y2=chirp(t,(f(i)-15/2)*carrier_freq_separation,switch_time,(f(i+1)-15/2)*carrier_freq_separation,'linear',90);
        y=y1-y2*1i;
        y =[y, y];
    end
    carrier_num=f(Nhop);
    sim('fsk2_v2');
    y=[y, reshape(signal_out,1,[])];
elseif strcmp(modulation_mode,'bpsk')
    sim('bpsk_v2');
end