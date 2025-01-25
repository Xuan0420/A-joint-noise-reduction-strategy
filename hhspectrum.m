function [A,f,tt] = hhspectrum1(x,t,l,aff)
% [A,f,tt] = HHSPECTRUM(imf,t,l,aff) computes the Hilbert-Huang spectrum
%
% inputs:
% 	- imf : matrix with one IMF per row
%   - t   : time instants
%   - l   : estimation parameter for instfreq
%   - aff : if 1, displays the computation evolution
%
% outputs:
%   - A   : amplitudes of IMF rows
%   - f   : instantaneous frequencies
%   - tt  : truncated time instants
%
% calls:
%   - hilbert  : computes the analytic signal
%   - instfreq : computes the instantaneous frequency
error(nargchk(1,4,nargin));
if nargin < 2
    t=1:size(x,2);
end
if nargin < 3
    l=1;
end
if nargin < 4
    aff = 0;
end
if min(size(x)) == 1
    if size(x,2) == 1
        x = x';
        if nargin < 2
            t = 1:size(x,2);
        end
    end
    Nmodes = 1;
else
    Nmodes = size(x,1);
end
lt=length(t);
tt=t((l+1):(lt-l));
for i=1:Nmodes
    an(i,:)=hilbert(x(i,:)')';
    f(i,:)=instfreq(an(i,:)',tt,l)';
    A=abs(an(:,l+1:end-l));
    if aff
        disprog(i,Nmodes,max(Nmodes,100))
    end
end
end

