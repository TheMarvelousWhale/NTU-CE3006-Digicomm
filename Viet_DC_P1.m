clear all; 
close all; 
clc;

%---------------Define Declaration---------------%
%Num of bits = 1024
nBits = 1024;

%Signal power = 1                   
Signal_Power = 1; 

%SNR_dB = 10log(Signal_Power/Noise_Power)
%Create matrix to store SNR
SNR_dB = 0:0.5:10;

%SNR = S/N = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

%Set run times
Total_Run = 20;

%Iterate through different SNR value
for i = 1 : length(SNR)
	Avg_Error = 0;
	for j = 1 : Total_Run
		%Generate random binary digits(0 or 1), INPUT SIGNAL
		Data = round(rand(1,nBits));
        
		%Convert binary digit to (-1 or +1)
		Signal = 2 .* Data - 1;
        
		%Generate equal number of noise samples  		
		Noise_Power = Signal_Power ./SNR(i);
        
        %The randn() function is for normal distribution
		Noise = sqrt(Noise_Power/2) .*randn(1,nBits);   
        
		%Output Signal
		Receive = Signal+Noise;                      

		Threshold = 0;
		Error = 0;
        
        %Fix the threshold value as 0
        %If received signal >= threshold value, threshold = 1
        %If received signal < threshold value, threshold = 0
		for k= 1 : nBits
			if (Receive(k)>= Threshold) && Data(k)==0||(Receive(k)<Threshold && Data(k)==1)
				Error = Error+1;
			end
        end
        
		%Calculate bit error rate during transmission
		Error = Error ./nBits;  
        
		%Calculate the average error for every runtime		
		Avg_Error = Error + Avg_Error;                   
	end
	Error_Rate(i) = Avg_Error / Total_Run;
end

%Calculate analytical Bit Error Rate
Pred_Rate=(1/2)*erfc(sqrt(SNR)); 

%Graph and Plot the result           
figure(1)
semilogy (SNR_dB,Error_Rate,'r*');
xlabel('Normalized SNR')
ylabel('Probability Error');

hold on
semilogy (SNR_dB,Pred_Rate,'m');
legend('simulation','theory');
axis([0 20 10^(-5) 1]);
hold off

%data generation
figure("position", [10,100,700,300]) 
plot(Signal);
title('Data_Generated')
%noise generation
figure("position", [10,500,700,300]) 
plot(Noise);
title('Noise_Generated')
%received data generation
figure("position", [800,100,700,300]) 
plot(Receive);
title('Received_Data')



