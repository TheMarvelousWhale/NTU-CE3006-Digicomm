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

%calculate signal length
SignalLength = fs*nBits/dataRate + 1;

%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = 0:1:20;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

%set run times
Total_Run = 10;

%define placeholder for error calculation
Error_RateOOK = zeros(length(SNR));
Error_RateBPSK = zeros(length(SNR));
Error_RateBFSK = zeros(length(SNR));

%FOR BFSK
deltaf = .5;
fh= fc + (fc*deltaf);
fl= fc - (fc*deltaf);
carh = cos(2*pi*fh*t)
carl = cos(2*pi*fl*t)

%BANDPASS FOR BFSK
%other value [0.1 0.2]
[b,a] = butter(6, [0.18 0.2], 'bandpass');

%for each SNR value
for i = 1 : length(SNR)
	Avg_ErrorOOK = 0;
    Avg_ErrorBPSK = 0;
    Avg_ErrorBFSK = 0;
    
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

        
        %%%%%%%%%%%%%%%%%%%BFSK%%%%%%%%%%%%%%%%%%%%%%
        %Generate BFSK Signal
        DataStreamBFSKOnes = DataStream .* Amp .* carh;
        DataStreamBFSKZeros = (DataStream - 1) .* Amp .* -carl;
        SignalBFSK = DataStreamBFSKOnes + DataStreamBFSKZeros;
        
        %BFSK Signal Power
        Signal_Power_BFSK = (norm(SignalBFSK)^2)/SignalLength;
        
        %Generage Noise BFSK
        Noise_Power_BFSK = Signal_Power_BFSK ./SNR(i);
        NoiseBFSK = sqrt(Noise_Power_BFSK/2) .*randn(1,SignalLength);
        
        %Receive BFSK Signal
        ReceiveBFSK = SignalBFSK+NoiseBFSK;
        
        %Demodulation
        %BandPass Filter
        FilteredBFSKZeros = filtfilt(b, a, ReceiveBFSK);
        FilteredBFSKOnes = filtfilt(b, a, ReceiveBFSK);
        
        %Envelope Detector
        [upperenvBFSKZeros, lowerenvBFSKZeros] = envelope(FilteredBFSKZeros, 'linear');
        [upperenvBFSKOnes, lowerenvBFSKOnes] = envelope(FilteredBFSKOnes, 'linear');
        
        %Summation
        FilteredBFSK = upperenvBFSKOnes + upperenvBFSKZeros;
        
        %Sample and Detect
        samplingPeriod = fs / dataRate;
        sampledBFSK = sample(FilteredBFSK, samplingPeriod, nBits);
        
        resultBFSK = decision_device(sampledBFSK,nBits,Amp/2);
        %%%%%%%%%%%%%%%END OF BFSK%%%%%%%%%%%%%%%
        
        
        %----- OOK -----%
        Signal_OOK = Carrier .* DataStream;
        
        %generate noise 
        Signal_Power_OOK = (norm(Signal_OOK)^2)/SignalLength;
		Noise_Power_OOK = Signal_Power_OOK ./SNR(i);
		NoiseOOK = sqrt(Noise_Power_OOK/2) .*randn(1,SignalLength);
		
        %transmission
		ReceiveOOK = Signal_OOK+NoiseOOK;
        %detection
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
        
        %Multiple and Low Pass Filter
        MultipliedBPSK = DividedBPSK .* ReceiveBPSK;
        OutputBPSK = filtfilt(b_low, a_low, MultipliedBPSK);
        
        %sample and decision device
        sampledBPSK = sample(OutputBPSK, samplingPeriod, nBits);
        resultBPSK = decision_device(sampledBPSK,nBits,0);           %-- bipolar -- threshold 0        
        
        
        %--Calculate Error--%
        ErrorOOK = 0;
        ErrorBPSK = 0;
        ErrorBFSK = 0;
        
        for k = 1: nBits - 1
            if(result_OOK(k) ~= Data(k))
                ErrorOOK = ErrorOOK + 1;
            end
            if(resultBPSK(k) ~= Data(k))
                ErrorBPSK = ErrorBPSK + 1;
            end
            if(resultBFSK(k) ~= Data(k))
                ErrorBFSK = ErrorBFSK + 1;
            end
        end
        Avg_ErrorOOK = ErrorOOK + Avg_ErrorOOK;
        Avg_ErrorBPSK = ErrorBPSK + Avg_ErrorBPSK;
        Avg_ErrorBFSK = ErrorBFSK + Avg_ErrorBFSK;
	end
	Error_RateOOK(i) = Avg_ErrorOOK / Total_Run;
    Error_RateBPSK(i) = Avg_ErrorBPSK / Total_Run;
    Error_RateBFSK(i) = Avg_ErrorBFSK / Total_Run;
end

%Error plot
figure(1);
semilogy (SNR_dB,Error_RateOOK,'k-*');
hold on
semilogy(SNR_dB, Error_RateBPSK, 'c-*');
hold on
semilogy(SNR_dB, Error_RateBFSK, 'y-+');
hold off
title('Error rate of OOK and BPSK and BFSK for different SNR');
legend('OOK', 'BPSK', 'BFSK');
ylabel('Pe');
xlabel('Eb/No')

%OOK plot
figure(2);
subplot(221);title('Transmitted Signal OOK');plot(Signal_OOK,'k');
subplot(222);title('Received Signal OOK');plot(ReceiveOOK, 'k')
subplot(223);title('Filtered Signal OOK');plot(FilteredOOK, 'k');
subplot(224);title('Captured Data');plot(sampledOOK);

%BPSK plot
figure(3)
subplot(221);title('Transmitted Signal BPSK');plot(Signal_BPSK,'k');
subplot(222);title('Received Signal BPSK');plot(ReceiveBPSK, 'k')
subplot(223);title('Filtered Signal BPSK');plot(FilteredBPSK, 'k');
subplot(224);title('Captured Data');plot(sampledBPSK);

%BFSK plot
figure(4)
subplot(221);title('Transmitted Signal BFSK');plot(SignalBFSK,'k');
subplot(222);title('Received Signal BFSK');plot(ReceiveBFSK, 'k')
subplot(223);title('Filtered Signal BFSK');plot(FilteredBPSK, 'k');
subplot(224);title('Captured Data');plot(sampledBFSK);

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

function [upperenv lowerenv] = envelope(sig, method)
if nargin == 1 
    method = 'linear';
end
upperind = find(diff(sign(diff(sig))) < 0) + 1;
lowerind = find(diff(sign(diff(sig))) > 0) + 1;
f = 1;
l = length(sig);
try
    upperind = [f upperind l];
    lowerind = [f lowerind l];
catch 
    upperind = [f; upperind; l];
    lowerind = [f; lowerind; l];
end
xi = f : l;
upperenv = interp1(upperind, sig(upperind), xi, method, 'extrap');
lowerenv = interp1(lowerind, sig(lowerind), xi, method, 'extrap');
end