function [Ug] = create_signal(Ag, Fs)

%Constants
Xj = .3; %vocal fold thickness
P0 = 1.14 * 10^-3; %density
Mu = 1.86 * 10^-4; %viscosity coefficient
Lj = 1.2; %vocal fold length
kc = 1.42;
Psub = 8*980.39;

Lg = zeros(size(Ag));
Rg = zeros(size(Ag));
Q = zeros(size(Ag));
Ug = zeros(size(Ag));

Q(1) = 0;
Ug(1) = 0;
Lg(1) = P0 * Xj / Ag(1);
Rg(1) = 4.0 * pi * Mu / Ag(1);

for u = 2:(size(Ag)(1))
    if Ag(u) ~= 0
        Lg(u) = P0 * Xj / Ag(u);
        Rg(u) = ((12 * Mu * Lj^2) / (Ag(u)^3)) + (kc * P0 * Ug(u-1) / (Ag(u)^2));
        Q(u) = Fs * Lg(u-1) * Ug(u-1) * Q(u-1);
        Ug(u) = (Psub + Q(u)) / (2 * Fs * Lg(u) + Rg(u));
    end
end

sound(diff(Ug), Fs);

t=(1/Fs):(1/Fs):D;
plot(1000*t(1:1000),Ug(1:1000));
xlabel('milliseconds');


end

