function [Agp] = glottal_area(Ap, F0, Fs, D)
% Ap: peak amplitude
% F0: fundamental frequency
% Fs: sampling rate
% D: duration of sample

Agp = zeros(Fs * D, 1);

T0 = 1.0 / F0;
t1 = .36 * T0;
t2 = .26 * T0;
a = pi / t1;
b = 1.0 / (1.0 - cos(pi * t2 / t1));
dt = 1.0 / Fs;

for n = 1:(Fs * D)
    t = mod((n * dt), T0);
    if t < t1
        Agp(n) = .5 * Ap * (1.0 - cos(a * t));
    elseif t < (t1 + t2)
        Agp(n) = Ap * (1.0 - b + b * cos(a * (t - t1)));
    else
        Agp(n) = .00001;
    end
end
Agp(Fs*D) = .00001;

%figure;
%t = dt:dt:D;
%plot(t, Agp);

end

