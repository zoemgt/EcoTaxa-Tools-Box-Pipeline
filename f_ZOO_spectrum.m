% f_ZOO_spectrum.m
% Computes a size spectrum from ZooScan / imaging data.
%
% Inputs:
%   sss  : log(biovolume) values for each organism [log(mm3)]
%   vol  : biovolume values for each organism [mm3]
%   smin : lower biovolume limit
%   smax : upper biovolume limit
%   kk   : logarithmic bin base (e.g. 2^(1/4))
%   opt  : 1 = abundance spectrum | 2 = biovolume spectrum
%
% Outputs:
%   X  : bin centres (log scale)
%   Y  : normalised spectrum (abundance or biovolume per bin width)
%   Z  : unnormalised counts
%   X1 : bin widths

function [X, Y, Z, X1] = f_ZOO_spectrum(sss, vol, smin, smax, kk, opt)

%% Build log-spaced bin edges
xvec = smin;
while xvec(end) < smax
    xvec(end+1, 1) = smin * kk^(numel(xvec));
end

X1   = diff(xvec);
xlog = log(xvec);
X    = (xlog(1:end-1) + xlog(2:end)) / 2;  % bin centres

%% Binned counts
if opt == 1
    y1 = ones(numel(sss), 1);   % abundance: count each organism once
else
    y1 = vol;                   % biovolume: sum biovolumes per bin
end

n = numel(xlog);
Y = zeros(n-1, 1);
for t = 1:n-1
    in_bin = sss >= xlog(t) & sss < xlog(t+1);
    Y(t)   = sum(y1(in_bin));
end

Z = Y;              % unnormalised
Y = Y ./ X1;        % normalise by bin width

end
