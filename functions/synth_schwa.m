%function signal = synth_schwa(Area, Length, Ap, F0, Fs, duration)
A = 1; % area function
X = 17; % length of section
Ap = 0.2;
F0 = 100;
Fs = 10000;
duration = 1;

rho = 1.14 * 10^-3;
c = 3.5 * 10^4;
mu = 1.86 * 10^-4;
bval = 1.4*10^3;
m = 1.5;
k = 3 * 10^5;
Psub = 8*980.39;
T = 1/Fs;
glt_len=1.2;
glt_thck=0.3;
kc = 1.42;

[Agp] = glottal_area(Ap, F0, Fs, duration);
Agp = max(Agp,0.001); % To avoid divison be zero later

N=length(Agp); % this is how many samples I will have

Ulips = zeros(N,1); % this will be the output

Ud = 0; 
Q1 = 0; Q2 = 0; %initial value; it will change
V1 = 0; V2 = 0;
Vc = 0; Qwl = 0; Qwc = 0;

U(1)=1;
U(2)=2;

for n=2:N
    
    % j=1 calculate
    
    Lg = (2/T) * rho * glt_thck / Agp(n);
    Rg = (12 * mu * glt_len^2 * glt_thck) / (Agp(n)^3) + ((1.38/2) * rho * abs(U(1)) / (Agp(n)^2));
    
    L = (1/T) * rho * X / A;
    R = 4 * pi * mu  / (X * A^2);
    C = (2/T) * X * A / (rho * c^2);
    Lw = 2 * m / (T * A * X * sqrt(pi));
    Rw = bval / (A * X * sqrt(pi));
    Cw = T/(2*sqrt(pi)) * k / (A * X);
    Yw = 1 / (Lw + Cw + Rw);
    b1 = 1 / (C + Yw);

    Ud = 0; % Because Area and Length are constant over time
    
    % j=1 assign
    
    H1 = -(Lg + L + Rg + R + b1);
    
    % j=2 calculate
    
    Grad = 9 * pi^2 * A / (128 * rho * c);
    Srad = T*(1.5 * pi * sqrt(pi * A)) / (8 * rho);
    b2 = 1 / (Srad + Grad);
    
    H2 = -(b2 + b1 + R + L);

    % refresh force constants
    
    V1 = Vc - Yw * (Qwl - Qwc);
    F1 = - b1 * (Ud - V1) - Psub - Q1; 
    F2 = b1 * (Ud - V1) + b2 * V2 - Q2;
    
    % Create linear system and solve
    H = [H1, b1; b2, H2];
    F = [F1; F2];
    U = inv(H) * F;
    Ulips(n) = U(2);
    
    % glottis
    Q1 = 2 * (Lg + L) * U(1) - Q1; 
    P1 = b1 * (U(1) - U(2) - Ud + V1);
    
    
    Vc = 2 * C * P1 - Vc;
    u3 = Yw * (P1 + Qwl - Qwc);
    Qwl = 2 * Lw * u3 - Qwl;
    Qwc = 2 * Cw * u3 + Qwc;
    
    % lips
    P2 = b2 * (U(2) + V2);
    Q2 = 2 * L * U(2) - Q2;
    V2 = -2 * Srad * P2 + V2;
    
end

signal = diff(Ulips);% (diff(U2) - mean(diff(U2)))/std(diff(U2));

plot(signal(1:500));
soundsc(signal, Fs);