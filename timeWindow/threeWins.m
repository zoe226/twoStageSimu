%{
    hanning
    hamming
    rectangle
%}
clc
clear
close all

win_len = 256;
% matlab
hann_m = hanning(win_len).';
hamm_m = hamming(win_len).';
rec_m = boxcar(win_len).';

% compute
hann_c = zeros(1,win_len);
hamm_c = zeros(1,win_len);
rec_c = ones(1,win_len);

for i = 1:win_len
    hann_c(i) = 0.5*(1-cos(2*pi*(i-1)/(win_len-1)));
end

for i = 1:win_len
    hamm_c(i) = 0.54 - 0.46*cos(2*pi*(i-1)/(win_len-1));
end

figure;
plot(hann_m);hold on;plot(hann_c);hold off;
title('hanning');legend('hannM','hannC');
figure;
plot(hamm_m);hold on;plot(hamm_c);hold off;
title('hamming');legend('hammM','hammC');
figure;
plot(rec_m);hold on;plot(rec_c);hold off;
title('rectangle');legend('recM','recC');