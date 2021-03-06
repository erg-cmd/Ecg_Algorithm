function Hd = FIR_Equiri
%FIR_EQUIRI Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 8.3 and the DSP System Toolbox 8.6.
% Generated on: 09-Jun-2015 02:50:03

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 360;  % Sampling Frequency

Fpass = 14;              % Passband Frequency
Fstop = 70;              % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.00001;          % Stopband Attenuation
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]
