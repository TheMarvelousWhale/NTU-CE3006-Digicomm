%--Admin stuff--%
clear all; close all; clc;


%define carrier frequency
fc = 10000; %10kHz
%16 times oversampled -> sample freq = 16 fc
fs = 16 * fc;

%define data rate of 1kbps
dataRate = 1000;
%define number of data bits
nBits = 1024;
%define sampling rate
samplingPeriod = fs / dataRate;

%define Amplitude
Amp = 5;
%define time steps
t = 0: 1/fs : nBits/dataRate;

%define 6th order LP butterworth filter with 0.2 normalized cutoff frequency
[b_low,a_low] = butter(6, 0.2);
%define 6th order HP butterworth filter with 0.2 normalized cutoff frequency
[b_high,a_high] = butter(6, 0.2, 'high');


%generate carrier frequency
Carrier = Amp .* cos(2*pi*fc*t);
Carrier_2 = Amp .* cos(2*pi*10*fc*t);
%calculate signal length
SignalLength = fs*nBits/dataRate + 1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = -20:1:20;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

%MODIFY THE VARIABLE BELOW TO CHOOSE AT WHICH SNR VALUE 
%TO PLOT SIGNAL,NOISE and RECEIVE
plot_SNR_dB = 15;



%set run times
Total_Run = 10;

%define placeholder for error calculation
Error_RateOOK = zeros(length(SNR));
Error_RateBPSK = zeros(length(SNR));
Error_RateFSK = zeros(length(SNR));


%for each SNR value
for i = 1 : length(SNR)
	Avg_ErrorOOK = 0;
    Avg_ErrorBPSK = 0;
    Avg_ErrorFSK = 0;
    
    
    %for each SNR value, average the error over %Total_Run times
	for j = 1 : Total_Run
        
        %-----Data generation-----%
        Data = round(rand(1,nBits));
        
        %fill the data stream
        DataStream = zeros(1, SignalLength);
        for k = 1: SignalLength - 1
            DataStream(k) = Data(ceil(k*dataRate/fs));
        end
        DataStream(SignalLength) = DataStream(SignalLength - 1);

        
        
        %----- OOK -----%
        Signal_OOK = Carrier .* DataStream;
        
        %generate noise 
        Signal_Power_OOK = (norm(Signal_OOK)^2)/SignalLength;  %Sum of squared signal amp over signal length
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
        
        
        
        %----- binary phase shift keying -----%
        DataStream_BPSK = DataStream .* 2 - 1;
        Signal_BPSK = Carrier .* DataStream_BPSK;
        
        %generate noise
        Signal_Power_BPSK = (norm(Signal_BPSK)^2)/SignalLength;
		Noise_Power_BPSK = Signal_Power_BPSK ./SNR(i);
		NoiseBPSK = sqrt(Noise_Power_BPSK/2) .*randn(1,SignalLength);
		
        %transmission
        ReceiveBPSK = Signal_BPSK+NoiseBPSK;
        %non-coherent detection -- square law
        SquaredBPSK = ReceiveBPSK .* ReceiveBPSK;
        %high pass filter 
        FilteredBPSK = filtfilt(b_high, a_high, SquaredBPSK);
    
        %frequency divider
        DividedBPSK = interp(FilteredBPSK, 2);
        DividedBPSK = DividedBPSK(1:length(FilteredBPSK));
        
        %Multiple and Band Pass Filter
        MultipliedBPSK = DividedBPSK .* ReceiveBPSK;
        OutputBPSK = filtfilt(b_low, a_low, MultipliedBPSK);
        
        %sample and decision device
        sampledBPSK = sample(OutputBPSK, samplingPeriod, nBits);
        resultBPSK = decision_device(sampledBPSK,nBits,0);           %-- bipolar -- threshold 0        
        
        
        %--FSK--%
        DataStreamFSK = -1.*DataStream + 1; %bit flip to multiply it with a diff carrier
        Signal_FSK = Carrier_2 .* DataStream + Carrier.*DataStreamFSK; 
        
        %generate noise
        Signal_Power_FSK = (norm(Signal_FSK)^2)/SignalLength;
		Noise_Power_FSK = Signal_Power_FSK ./SNR(i);
		NoiseFSK = sqrt(Noise_Power_FSK/2) .*randn(1,SignalLength);
        
        %transmission
		ReceiveFSK = Signal_FSK+NoiseFSK;
        %detection -- bandpass filters
        ReceiveFSK_LOW = filtfilt(b_low,a_low,ReceiveFSK);
        ReceiveFSK_HIGH = filtfilt(b_high,a_high,ReceiveFSK);
        %Envelope 
        SquaredFSK_LOW = ReceiveFSK_LOW.*ReceiveFSK_LOW;
        SquaredFSK_HIGH = ReceiveFSK_HIGH.*ReceiveFSK_HIGH;
        SquaredFSK = SquaredFSK_HIGH - SquaredFSK_LOW;
        %sample and decision device
        SampledFSK = sample(SquaredFSK, samplingPeriod, nBits);
        resultFSK = decision_device(SampledFSK,nBits, 0);  
        
        

        %--Calculate Error--%
        ErrorOOK = 0;
        ErrorBPSK = 0;
        ErrorFSK = 0;
        
        for k = 1: nBits - 1
            if(result_OOK(k) ~= Data(k))
                ErrorOOK = ErrorOOK + 1;
            end
            if(resultBPSK(k) ~= Data(k))
                ErrorBPSK = ErrorBPSK + 1;
            end
            if(resultFSK(k) ~= Data(k))
                ErrorFSK = ErrorBPSK + 1;
            end
        end
        Avg_ErrorOOK = ErrorOOK + Avg_ErrorOOK;
        Avg_ErrorBPSK = ErrorBPSK + Avg_ErrorBPSK;
        Avg_ErrorFSK = ErrorFSK + Avg_ErrorFSK;
    end
    if (plot_SNR_dB == SNR_dB(i))
            plot_signal = Data;
            plot_mod_OOK = Signal_OOK;
            plot_receive_OOK = ReceiveOOK;
            plot_demod_OOK = FilteredOOK;
            plot_decoded_OOK = result_OOK;
            plot_mod_BPSK = Signal_BPSK;
            plot_receive_BPSK = ReceiveBPSK;
            plot_demod_BPSK = FilteredBPSK;
            plot_decoded_BPSK = resultBPSK;
            plot_mod_FSK = Signal_FSK;
            plot_receive_FSK = ReceiveFSK;
            plot_demod_FSK = SquaredFSK;
            plot_decoded_FSK = resultFSK;
    end
    
	Error_RateOOK(i) = (Avg_ErrorOOK / Total_Run)/nBits;
    Error_RateBPSK(i) = (Avg_ErrorBPSK / Total_Run)/nBits;
    Error_RateFSK(i) = (Avg_ErrorFSK / Total_Run)/nBits;
end


%Error plot
figure(1);
semilogy (SNR_dB,Error_RateOOK,'k-*');
hold on
semilogy(SNR_dB, Error_RateBPSK, 'c-*');
hold off
hold on
semilogy(SNR_dB, Error_RateFSK, 'r-*');
hold off
title('Error rate of OOK and BPSK for different SNR');
legend('OOK', 'BPSK','FSK');
ylabel('Pe');
xlabel('Eb/No')


%OOK plot
figure(2);
subplot(511);title('Generated Data');plot(plot_signal);
subplot(512);title('Modulated OOK');plot(plot_mod_OOK,'k');
subplot(513);title('Received Signal OOK');plot(plot_receive_OOK, 'k')
subplot(514);title('Demodulated OOK');plot(plot_demod_OOK, 'k');
subplot(515);title('Decoded Data');plot(plot_decoded_OOK);

%BPSK plot
figure(3)
subplot(511);title('Generated Data');plot(plot_signal);
subplot(512);title('Modulated BPSK');plot( plot_mod_BPSK,'k');
subplot(513);title('Received Signal BPSK');plot(plot_receive_BPSK, 'k')
subplot(514);title('Demodulated BPSK');plot(plot_demod_BPSK, 'k');
subplot(515);title('Decoded Data');plot(plot_decoded_BPSK);

%FSK
figure(4)
subplot(511);title('Generated Data');plot(plot_signal);
subplot(512);title('Modulated FSK');plot( plot_mod_FSK,'k');
subplot(513);title('Received Signal FSK');plot(plot_receive_FSK, 'k')
subplot(514);title('Demodulated FSK');plot(plot_demod_FSK, 'k');
subplot(515);title('Decoded Data');plot(plot_decoded_FSK);

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