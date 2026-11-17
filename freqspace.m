% freqspace, version May 2011
% Software developed by
% LabMath-Indonesia & University of Twente
function [k,N] = freqspace(x,N,kind)
% FREQSPACE  Wavenumber (or frequency) vector for a uniformly spaced grid x.
%   k = freqspace(x)         FFT-ordered wavenumbers for N = length(x) points
%   k = freqspace(x,N,kind)  kind: 'default' (FFT order), 'unshifted'
%                            (monotonic), 'single sided'
%   NOTE: this shadows MATLAB's built-in freqspace (Image Processing Toolbox).

if nargin == 1
    N    = length(x);
    kind = 'default';
elseif nargin == 2
    kind = 'default';
end

dx   = x(2)-x(1);
L    = x(end)-x(1);
dk   = 2*pi/L;
fmin = -0.5*2*pi/dx;
if rem(N,2) == 0
    fmax = -fmin-dk;
elseif rem(N,2) == 1
    fmax = -fmin;
end

switch kind
    case 'default'
        k = ifftshift(linspace(fmin,fmax,N)');
    case 'single sided'
        k = linspace(0,-fmin,(floor(N/2)+1)*2*pi);   % factor 2pi
    case 'unshifted'
        k = linspace(fmin,fmax,N)';
end
end
