close all
clear sound
N = 100000000; % number of samples

% note the term 'sum' is interchangable with 'left' ect. depending on what
% information is presented in the imported files

file_input_sum = fopen('left.dat','r');
input_sum = fread(file_input_sum,N,'float'); % read for n files
fclose(file_input_sum);
file_input_diff = fopen('right.dat','r');
input_diff = fread(file_input_diff,N,'float'); % read for n files
fclose(file_input_diff);

fs = 48e3; % sample rate in Hz   CHANGED

t1 = [1:length(input_sum)]'/fs;
fl = length(t1);
fl = 2^ceil(log2(fl));
f = (-fl/2:fl/2-1)/(fl/fs);
SUMFFT = fftshift(fft(input_sum,fl));

t = [1:length(input_diff)]'/fs;
fl = length(t);
fl = 2^ceil(log2(fl));
f2 = (-fl/2:fl/2-1)/(fl/fs);
DIFFFFT = fftshift(fft(input_diff,fl));


figure(1);
subplot(2,2,1);
plot(input_sum);
title('Time Response of LEFT');
xlabel('samples');
ylabel('LEFT');
grid;

subplot(2,2,2);
plot(input_diff);
title('Time Response of RIGHT');
grid;
xlabel('samples');
ylabel('RIGHT');

subplot(2,2,3);
plot(f,abs(SUMFFT));
title('Freq Responce of LEFT');
xlabel('f(Hz)');
ylabel('Amplitude');
grid;

subplot(2,2,4);
plot(f,abs(DIFFFFT));
title('Freq Responce of RIGHT');
xlabel('f(Hz)');
ylabel('Amplitude');
grid;

disp('Finished')

soundsc([input_sum , input_diff],48000)

% input = input/(max(abs(input))); % normalise output to amplitude 1
% audiowrite('output.wav',input,48000);