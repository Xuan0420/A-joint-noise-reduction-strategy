function psd = psd_entropy(modes, fs)
    % Calculate the power spectral entropy of each mode component and the original noise-free signal
    for i = 1:size(modes, 1)
        % Calculate the power spectral density of the mode component
        [pxx, ~] = pwelch(modes(i, :));
        % Calculate the power spectral entropy of the mode component
        psd(i) = entropy(pxx);
    end
end
