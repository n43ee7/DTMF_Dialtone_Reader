%% Nabeel Nayyar | DFT Project | CSE 3313 Introduction to Signal Processing
%
% Project: Design a filter to identify frequencies of DTMF Dial tones
%           
%
% Notes: The following DTMF Codebase was developed using .wav samples 
%        generated from an online resourse available here
%        (https://www.audiocheck.net/audiocheck_dtmf.php)
%        
%       
%% Program House Keeping
close all;
clear;
clearvars;
clc;
warning off; % Believe me it works

% Local running path
oldpth = cd;

%% Runtime Parameters
% Adjust ONLY if chaning the type of input signal
GRAIN_DURATION = 50*(10^-3);    % Signal Grain chop duration [Current: 50ms]
BASE_SIGNAL_GAIN = 5;           % Signal Power Selection threshold [Current: 5dB]
DTMF_ADJ_F = 10;                % DTMF Frequency Comparison tolerance [Current: 10 Hz]

%% User Input Collection
% We begin by reading the input wave file to the system
[file, path] = uigetfile('*.wav');
if isequal(file, 0)
   disp('User: Cancel');
   return; % Exit
else
   disp(['User Select: ', fullfile(path,file)]);
end

cd(path);
[data, SampleRate] = audioread(file);
soundsc(data, SampleRate)
s_size = length(data);

%% Input Sample Detailing
% Dynamically adjusting a time matrix for input signal
time = (0:1:length(data)-1)*(1/SampleRate);

% Plotting the original Input Signal
figure
subplot(211);
plot(time, data)
xlabel('T [s]')
ylabel('Amp [V]')
title('Input Signal Waveform')

% Signal Spectrum Analysis
scale = length(data)+1;
% Linearly spaced vector of frequency response between the points(0-30400)/(30399)
f_resp = linspace(0, SampleRate, scale);
f_resp(end)=[]; % Trim last value to match data[NxM] and f_re[MxN]

% TESTING Taking DFT of the input signal with calculated scale
dft_data = fft(data, scale-1);

subplot(212)
plot(f_resp/SampleRate, abs(dft_data))
xlabel('f [kHz]')
ylabel('PSD')
title('Input Signal DFT')

%% Granular synthesis of Input sample
% Identifying individual sample grains and checking any activity in signal

% Grains = chop duration / sample period 
s_period = 1/SampleRate;
G_rate = GRAIN_DURATION/s_period;

% Chopping the signal in grain size for later frequency analysis
% Round each element of X to the nearest integer greater than or equal to that element 
G_dur = ceil(s_size/G_rate);
G_dur = G_dur * G_rate - s_size;

% Transposed data padded with zeros to for adjusted frequency
pre_dft_data = transpose(reshape([transpose(data) zeros(1, G_dur)]', G_rate, []));
 
% Taking metrics for vector allocation
s_data_lim = size(pre_dft_data, 1);  % for a [MxN] Block >> M size


%% Application of DFT for filtering Results from Grains
% Calculated frequency matrix of the sample signal
f_mat = 0:(1 / s_size):(s_size - 1)*(1 / s_size);

% We start processing for the first value of the pre dft data 
% Applying DFT on transposed and adjusted audio data
% and adjusting the  Magnitudes by taking absolute of DFT
abs_dft_data_N = abs(fft(pre_dft_data(1,:), s_size));
abs_dft_data_N = abs_dft_data_N(1:(s_size + 1) / 2);

% Find peaks in the dft signal. This will give us the frequency of dialed
% samples in the signal input
[pks, locs] = findpeaks(abs_dft_data_N, 'SortStr', 'descend');

% From the frequency matrix of the signal we find the first peak 
Dial_Notes(1,1) = f_mat(locs(1));

% Setting a band pass filter for the grained signal for sorting frequencies
Temp_chebT = Dial_Notes(1,1);
Temp_chebT = [0.555 * Temp_chebT 0.9323 * Temp_chebT];
[b,a] = cheby1(7, 1, Temp_chebT, 'stop');

abs_dft_data_N = abs(fft(filter(b, a, pre_dft_data(1,:)), s_size));

[pks, locs] = findpeaks(abs_dft_data_N(1:(s_size + 1) / 2), 'SortStr', 'descend');
Dial_Notes(1, 2) = f_mat(locs(1));

% So this means we can now iterate till the limit of the data and repeat
% for all samples in the signal dft_pre matrix to find signal activation
% and store them in a Dialed Note matrix
for i = 2 : s_data_lim

    % We need to calculate the previous signals power to prevent alising
    G_pwr(i - 1) = sum(pre_dft_data(i - 1, :).^2) / G_rate;
    G_pwr_db(i - 1) = (10*log(G_pwr(i-1)/10.^-3));
       
    % As did above we calculate the power of the signal
    G_pwr(i) = sum(pre_dft_data(i, :).^2) / G_rate;
    G_pwr_db(i) = (10*log(G_pwr(i)/10.^-3));

    % Comparison of active signal gain to see if it is accepted and to
    % avoid signal aliasing of duplicate frequencies
    if G_pwr_db(i - 1) < BASE_SIGNAL_GAIN 
        if G_pwr_db(i) > BASE_SIGNAL_GAIN 
            % Magnitude of the Data after DFT
            abs_dft_data_N = abs(fft(pre_dft_data(i,:), s_size));
            % Finding Peaks and locating the frequency
            [pks, locs] = findpeaks(abs_dft_data_N(1:(s_size+1)/2),'SortStr','descend');
            % Adding the location of the found frequency in the Dialed Notes
            Dial_Notes(i,1) = f_mat(locs(1));
            Temp_chebT = Dial_Notes(i,1);

            % We do the same for finding the high-frequency band of Dialed Notes
            Temp_chebT = [0.555 * Temp_chebT 0.9323 * Temp_chebT];
            % Band-pass again 
            [b,a] = cheby1(7, 1, Temp_chebT, 'stop');

            % To calculate the Magnitude of the DFT of the adjusted data
            % for high band frequencies
            abs_dft_data_N = abs(fft(filter(b, a, pre_dft_data(i, :)),s_size)); 
            [pks, locs] = findpeaks(abs_dft_data_N(1:(s_size + 1) / 2), 'SortStr', 'descend');
            Dial_Notes(i,2) = f_mat(locs(1)); 
        end
    end
    
end

% the first and second column were sorted to make the first column to be
% the low frequency and second column the high frequency
Dial_Notes = reshape(nonzeros(Dial_Notes),[],2)*SampleRate;
Dial_Notes = transpose(sort(Dial_Notes.','ascend'));


%% Reading Data of the Dialed Numbers in Input Signal

% Showing Power spectrum density of the calculated signal
figure
N = numel(data);
t = (0:N-1)/SampleRate; 
subplot(2,1,1)
pspectrum(data,SampleRate,'Leakage',1,'FrequencyLimits',[650, 1500])

% Showing the dialed Tones and their low frequency and high frequency
subplot(212)
stem(1:length(Dial_Notes), Dial_Notes(:,1))
hold on;
stem(1:length(Dial_Notes), Dial_Notes(:,2))
xlabel('Tone Number #')
ylabel('F [Hz]')
title('Dailed Tone Frequencies')
legend('Low Frequency Response', 'High Frequency Response','Location','NorthOutside','Orientation','horizontal','Box','off')


%% Check Found DTMF Frequencies and the respective dialed Number
for i = 1:length(Dial_Notes)
    Dialed_Number(i) = string(DTMFAlloc(Dial_Notes, i, DTMF_ADJ_F));
end

% Finally Merging the string to one
Dialed_Number = join(Dialed_Number);


%% Results and Termination Routine 
fprintf("Dialed Number: %s\n\n", Dialed_Number);
cd(oldpth); % Adjust PWD for robust performance

%% DTMF Localizing Function
% Provided in the Project notes we are given a table of DTMF Frequencies we
% can use the input [2xN] Matrix where [1xN] is Low Band and [2xN] is High
% Band. The DTMF Assignment has custom tolerance for an input signal due to
% band pass attenuation or something else. Set accordingly

function DTMF_RES = DTMFAlloc(x , i, DTMF_ADJ_F)
    % >1	L:697	H:1209
    if x(i,1) >= 697-DTMF_ADJ_F && x(i,1) <= 697+DTMF_ADJ_F && x(i,2) >= 1209-DTMF_ADJ_F && x(i,2) <= 1209+DTMF_ADJ_F
        DTMF_RES = '1';
    % >2	L:697	H:1336
    elseif  x(i,1) >= 697-DTMF_ADJ_F && x(i,1) <= 697+DTMF_ADJ_F && x(i,2) >= 1336-DTMF_ADJ_F && x(i,2) <= 1336+DTMF_ADJ_F
        DTMF_RES = '2';
    % >3	L:697	H:1477
    elseif x(i,1) >= 697-DTMF_ADJ_F && x(i,1) <= 697+DTMF_ADJ_F && x(i,2) >= 1477-DTMF_ADJ_F && x(i,2) <= 1477+DTMF_ADJ_F
        DTMF_RES = '3';
    % >4	L:770	H:1209
    elseif x(i,1) >= 770-DTMF_ADJ_F && x(i,1) <= 770+DTMF_ADJ_F && x(i,2) >= 1209-DTMF_ADJ_F && x(i,2) <= 1209+DTMF_ADJ_F
        DTMF_RES = '4';
    % >5	L:770	H:1336
    elseif  x(i,1) >= 770-DTMF_ADJ_F && x(i,1) <= 770+DTMF_ADJ_F  && x(i,2) >= 1336-DTMF_ADJ_F && x(i,2) <= 1336+DTMF_ADJ_F
        DTMF_RES = '5';
    % >6	L:770	H:1477
    elseif x(i,1) >= 770-DTMF_ADJ_F && x(i,1) <= 770+DTMF_ADJ_F  && x(i,2) >= 1477-DTMF_ADJ_F && x(i,2) <= 1477+DTMF_ADJ_F
        DTMF_RES = '6';
    % >7	L:852	H:1209
    elseif x(i,1) >= 852-DTMF_ADJ_F && x(i,1) <= 852+DTMF_ADJ_F && x(i,2) >= 1209-DTMF_ADJ_F && x(i,2) <= 1209+DTMF_ADJ_F
        DTMF_RES = '7';
    % >8	L:852	H:1336
    elseif  x(i,1) >= 852-DTMF_ADJ_F && x(i,1) <= 852+DTMF_ADJ_F  && x(i,2) >= 1336-DTMF_ADJ_F && x(i,2) <= 1336+DTMF_ADJ_F
        DTMF_RES = '8';
    % >9	L:852	H:1477
    elseif x(i,1) >= 852-DTMF_ADJ_F && x(i,1) <= 852+DTMF_ADJ_F  && x(i,2) >= 1477-DTMF_ADJ_F && x(i,2) <= 1477+DTMF_ADJ_F
        DTMF_RES = '9';
    % >0	L:941	H:1336
    elseif  x(i,1) >= 941-DTMF_ADJ_F && x(i,1) <= 941+DTMF_ADJ_F && x(i,2) >= 1336-DTMF_ADJ_F && x(i,2) <= 1336+DTMF_ADJ_F
        DTMF_RES = '0';
    % >*    L:941	H:1209
    elseif x(i,1) >= 941-DTMF_ADJ_F && x(i,1) <= 941+DTMF_ADJ_F && x(i,2) >= 1209-DTMF_ADJ_F && x(i,2) <= 1209+DTMF_ADJ_F
        DTMF_RES = '*';
    % >#	L:941	H:1477
    elseif x(i,1) >= 941-DTMF_ADJ_F && x(i,1) <= 941+DTMF_ADJ_F  && x(i,2) >= 1477-DTMF_ADJ_F && x(i,2) <= 1477+DTMF_ADJ_F
        DTMF_RES = '#';    
    % > CANT BE READ
    else
        DTMF_RES = 'X';
    end
end