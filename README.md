# A-joint-noise-reduction-strategy


## Installation
Please copy the whole package under MATLAB platform. The key packages include signal processing toolbox, optimization toolbox, wavelet toolbox and statistics and machine learning toolbox.


## Usage
"SMAVMDIWTD.m" is the main function of the joint noise reduction strategy. For the purpose of illustration, people can run the main function directly. 
"EMDMain.m" is the function of empirical mode decomposition. It is used to process raw signals using empirical mode decomposition technique.
"VMDMain.m" is the function of variational mode decomposition. It is used to process raw signals using variational mode decomposition technique.
"VMD.m" is the function of variational mode decomposition. It is used for VMD processing of signals after optimization using SMA.
"BaoluoEntropy.m" is the envelope entropy function. It is used to calculate the envelope entropy of the signal.
"EvaMetrix.m" is the function that measures the effectiveness of noise reduction, including signal-to-noise ratio, peak signal-to-noise ratio, and root-mean-square error.
"hardThresholdDenoise.m" is the hard threshold function in traditional wavelet denoising. It can be used for hard threshold processing of signals.
"softThresholdDenoise.m" is the soft threshold function in traditional wavelet denoising. It can be used for soft threshold processing of signals.
"hhspectrum.m" is the function that computes computes the Hilbert-Huang Spectrum of a signal. The Hilbert-Huang Transform is a time-frequency analysis method that is often used for analyzing non-linear and non-stationary signals.
"hua_fft.m" defines a function hua_fft that performs the Fast Fourier Transform on a given signal y and provides options for plotting either the magnitude spectrum or the power spectrum of the signal.
"improved_threshold2.m" is the improved threshold function. It involves two parameters, alpha and beta, which modulate the behavior of the soft thresholding.
"improvedWaveletDenoise" is a signal denoising function based on wavelet transform and an improved soft thresholding technique.
"initialization.m" is a function that initializes a population of solutions for optimization problems.
"io.m" computes the index of orthogonality for a signal that has undergone empirical mode decomposition.
"objfun.m" is used to compute the minimum value of the local envelope entropy of the signal after performing a variational mode decomposition on the given noisy signal.
"psd_entropy.m" is a power spectrum entropy function, which is used to calculate the power spectrum entropy of each mode component and the original noiseless signal.
"SNR_singlech.m" is a function for calculating the signal-to-noise ratio of a noisy signal.
"PSNR_singlech.m" is a function of calculating the peak signal-to-noise ratio of a noisy signal.
"SMA.m" is the function of slime mold algorithm. The function optimizes a given objective function fobj by evolving a population of slime mold agents over several iterations. 
"threshold.m" is the main function of the wavelet threshold function. It is used to plot traditional soft and hard threshold functions and improved threshold functions.
"toimage.m" defines a function toimage that transforms a 2D matrix (A) into an image-like representation using the frequency information provided in f. It is designed for cases where the matrix (A) contains time-series data, and the goal is to reshape it into a frequency-time representation for visualization or further analysis.


## Remark
See the script file for more detailed comments.
