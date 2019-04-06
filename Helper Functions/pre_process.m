function y = pre_process(data)

global f fs nfft win;

noise_band = find(f>90,1):find(f<100,1,'last');                                 %NOISE?
temp_pow = pwelch(data,hanning(win),[],nfft,fs); % power spectral density estimate???
noise_lev = mean(temp_pow(noise_band));

data_smooth = smooth(data,10);
data_base = smooth(data_smooth,2000);
data_base = smooth(data_base,2000);
data_s = data_smooth - data_base;

y = data_s./sqrt(noise_lev);


