function y=signal_generate(hopping_rate, switch_time, carrier_freq_separation, modulation_mode, symbol_rate, SNR, signal_length, fs)
% Frequency-Hopping Signals Generater
% hopping_rate:                  跳频速率(hop/s)
% switch_time:                   切换时间
% carrier_freq_separation:       载波频率间隔(Hz)
% modulation_mode:               调制模式
% symbol_rate:                   符号速率
% SNR:                           信噪比
% signal_length:                 信号时间长度(s)
% fs:                            复信号采样率(Hz)

Nhop=fix(hopping_rate*signal_length);   %%仿真的跳数
hopping_cycle=1/hopping_rate; %跳频周期
symbol_cycle=1/symbol_rate; %符号周期
sample_per_symbol=fix(fs/symbol_rate);             %每一个符号的采样点数
sample_num_per_hop=fix(fs/hopping_rate);           %每一跳的采样点数
symbol_num_per_hop=fix(symbol_rate/hopping_rate);  %每一跳的符号数
sample_num_switch_time=fix(fs*switch_time);

f=zeros(Nhop,1);
for i=1:Nhop
    f(i)=randi(15,1,1);
end
y=[];
if strcmp(modulation_mode,'none')
    
elseif strcmp(modulation_mode,'2fsk')
    freq_separation=0.5*10^6; %2fsk调制频率间隔
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