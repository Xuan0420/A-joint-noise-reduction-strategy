function [cc, y_f] = hua_fft(y, fs, style, varargin)
% When style = 1, plot amplitude spectrum; when style = 2, plot power spectrum;
% when style = other values, plot both amplitude and power spectra.
% When style = 1, two additional optional parameters can be input.
% The optional parameters are used to control the frequency range to be displayed.
% The first parameter is the starting point of the frequency range.
% The second parameter is the ending point of the frequency range.
% Other styles do not have optional input parameters, if the input position is incorrect.

nfft = 2^nextpow2(length(y)); % Find the largest power of 2 greater than the length of y (automatically computes the best FFT step size nfft)
%nfft = 1024; % Manually set the FFT step size nfft
y = y - mean(y); % Remove DC component
y_ft = fft(y, nfft); % Perform DFT on y signal, obtaining the amplitude distribution in frequency
y_p = y_ft .* conj(y_ft) / nfft; % conj() function computes the conjugate complex of y, the conjugate of a real number is itself.
y_f = fs * (0:nfft/2-1) / nfft; % Frequency sequence corresponding to the transformed signal

if style == 1
    if nargin == 3
        cc = 2 * abs(y_ft(1:nfft/2)) / length(y);
        % ylabel('Amplitude'); xlabel('Frequency'); title('Signal Amplitude Spectrum');
        % plot(y_f, abs(y_ft(1:nfft/2))); % Method to plot FFT from the forum
    else
        f1 = varargin{1};
        fn = varargin{2};
        ni = round(f1 * nfft / fs + 1);
        na = round(fn * nfft / fs + 1);
        hold on
        plot(y_f(ni:na), abs(y_ft(ni:na) * 2 / nfft), 'k');
    end
elseif style == 2
    plot(y_f, y_p(1:nfft/2), 'k');
    % ylabel('Power Spectral Density'); xlabel('Frequency'); title('Signal Power Spectrum');
else
    subplot(211); plot(y_f, 2 * abs(y_ft(1:nfft/2)) / length(y), 'k');
    ylabel('Amplitude'); xlabel('Frequency'); title('Signal Amplitude Spectrum');
    subplot(212); plot(y_f, y_p(1:nfft/2), 'k');
    ylabel('Power Spectral Density'); xlabel('Frequency'); title('Signal Power Spectrum');
end
end
