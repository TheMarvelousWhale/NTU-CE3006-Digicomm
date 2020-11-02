%--Admin stuff--%
clear all; close all; clc;

fc = 10000; 
fs = 16 * fc;
dataRate = 1000;
nBits = 32;
samplingPeriod = fs / dataRate;   %for sampling later
Amp = 5;
t = 0: 1/fs : nBits/dataRate;    %steps is 1/fs ok; nbits/dataRate = total transmission duration
[b_low,a_low] = butter(6, 0.2);
[b_high,a_high] = butter(6, 0.2, 'high');
Carrier = Amp .* cos(2*pi*fc*t);
SignalLength = fs*nBits/dataRate + 1;                
SNR_dB = -20:1:20;
SNR = (10.^(SNR_dB/10));
plot_SNR_dB = 15;



%set run times
Total_Run = 10;

%define placeholder for error calculation
Error_RateOOK = zeros(length(SNR));
Error_RateBPSK = zeros(length(SNR));





for i = 1 : length(SNR)
	Avg_ErrorOOK = 0;
    Avg_ErrorBPSK = 0;
	for j = 1 : Total_Run
        Data = round(rand(1,nBits));
        DataStream = zeros(1, SignalLength);
        for k = 1: SignalLength - 1
            DataStream(k) = Data(ceil(k*dataRate/fs));
        end
        DataStream(SignalLength) = DataStream(SignalLength - 1);

        
        
        %----- OOK -----%
        Signal_OOK = Carrier .* DataStream;
        
        %generate noise 
        Signal_Power_OOK = (norm(Signal_OOK)^2)/SignalLength;
		Noise_Power_OOK = Signal_Power_OOK ./SNR(i);
		NoiseOOK = sqrt(Noise_Power_OOK/2) .*randn(1,SignalLength);
		
        %transmission
		ReceiveOOK = Signal_OOK+NoiseOOK;
        %detection -- square law device
        SquaredOOK = ReceiveOOK .* ReceiveOOK;
        %low pass filter
        FilteredOOK = filtfilt(b_low, a_low, SquaredOOK);
         
        %sample and decision device
        sampledOOK = sample(FilteredOOK, samplingPeriod, nBits);
        result_OOK = decision_device(sampledOOK,nBits, Amp/2);  %--OOK threshold is 0.5*(A+0)
        
        
        %--Calculate Error--%
        ErrorOOK = 0;
        ErrorBPSK = 0;
        
        for k = 1: nBits - 1
            if(result_OOK(k) ~= Data(k))
                ErrorOOK = ErrorOOK + 1;
            end
        end
        Avg_ErrorOOK = ErrorOOK + Avg_ErrorOOK;
        Avg_ErrorBPSK = ErrorBPSK + Avg_ErrorBPSK;

    end
    if (plot_SNR_dB == SNR_dB(i))
            plot_signal = Data;
            plot_mod_OOK = Signal_OOK;
            plot_receive_OOK = ReceiveOOK;
            plot_1 = SquaredOOK;
            plot_2 = FilteredOOK;
            plot_demod_OOK = FilteredOOK;
            plot_decoded_OOK = result_OOK;

    end
    
	Error_RateOOK(i) = (Avg_ErrorOOK / Total_Run)/nBits;
    Error_RateBPSK(i) = (Avg_ErrorBPSK / Total_Run)/nBits;
end





%OOK plot
figure(1);
subplot(611);title('Received Signal OOK');plot(plot_receive_OOK, 'k');
subplot(612);title('Squared OOK');plot(plot_1,'k');
subplot(613);title('Filtered OOK');plot(plot_2,'k');
subplot(614);title('Demodulated OOK');plot(plot_demod_OOK, 'k');
subplot(615);title('Decoded Data');plot(plot_decoded_OOK);
subplot(616);title('Data');plot(plot_signal);



%%--HELPER FUNCTION--%%
function sampled = sample(x,sampling_period,num_bit)
    sampled = zeros(1, num_bit);
    for n = 1: num_bit
        sampled(n) = x((2 * n - 1) * sampling_period / 2);
    end
end


%This function simulates the decision device
function binary_out = decision_device(sampled,num_bit,threshold)
    binary_out = zeros(1,num_bit);
    for n = 1:num_bit
        if(sampled(n) > threshold)
            binary_out(n) = 1;
        else 
            binary_out(n) = 0;
        end
    end
end