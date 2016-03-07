function signal = synth_static_vt_shape(A, X)
%(Area, Length, Ap, F0, Fs, duration)
M = length(A);
%A = zeros(M, 1); % area function; This should become an array with length M
%X = zeros(M, 1); % length of section; Also an array with length M
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
Q(1) = 0; Q(2) = 0; %initial value; it will change
V(1) = 0; V(2) = 0;
Vc(1) = 0; Qwl(1) = 0; Qwc(1) = 0;
[L, R, C, Lw, Rw, Cw, Yw, b, Ud, H, V, F, U, Q, P, Vc, u3, Qwl, Qwc] = deal(zeros(M+1,1));
U(1)=1;
U(2)=2;

for n=2:N
    
    % j=1 calculate
    
    Lg = (2/T) * rho * glt_thck / Agp(n);
    Rg = (12 * mu * glt_len^2 * glt_thck) / (Agp(n)^3) + ((1.38/2) * rho * abs(U(1)) / (Agp(n)^2));
    
    % This will loop over j, for j =1...M : all variables should become arrays, for
    % example the first one
    % L(j) = (1/T) * rho * X(j)/A(j);
    % Basically (almost) all variables that appera below will become
    % functions of j
for j = 1:M   
    L(j) = (1/T) * rho * X(j) / A(j);
    R(j) = 4 * pi * mu  / (X(j) * (A(j))^2);
    C(j) = (2/T) * X(j) * A(j) / (rho * c^2);
    Lw(j) = 2 * m / (T * A(j) * X(j) * sqrt(pi));
    Rw(j) = bval / (A(j) * X(j) * sqrt(pi));
    Cw(j) = T/(2*sqrt(pi)) * k / (A(j) * X(j));
    Yw(j) = 1 / (Lw(j) + Cw(j) + Rw(j));
    b(j) = 1 / (C(j) + Yw(j));      % b1->b(1)
end
    %Ud = 0; % Because Area and Length are constant over time
    
    % j=1 assign
    
    H(1) = -(Lg + L(1) + Rg + R(1) + b(1));
    % This will be for j = 1, i.e.
    % H(1) = -(Lg + L(1) + Rg + R(1) + b(1)
    % for j=2...M we will have
    % H(j) = - (L(j-1) + L(j) + R(j-1) + R(j) + b(j-1) + b(j))
for j = 2:M
    H(j) = - (L(j-1) + L(j) + R(j-1) + R(j) + b(j-1) + b(j));
end
    
    % j=2 calculate
    
    % Here A will become A(M)
    Grad = 9 * pi^2 * A(M) / (128 * rho * c);
    Srad = T*(1.5 * pi * sqrt(pi * A(M))) / (8 * rho);
    b(M+1) = 1 / (Srad + Grad); % This will be b(M+1)
    
    H(M+1) = -(b(M+1) + b(M) + R(M) + L(M)); % H(M+1) = -(b(M+1) + b(M) + R(M) + L(M))
    % Basically, wherever we had 2 in what follows will become M+1

    % refresh force constants
for j = 1:M  
    V(j) = Vc(j) - Yw(j) * (Qwl(j) - Qwc(j)); % V(1), same equation for j=1...M (but all vaules will be functions of j)
end

    F(1) = - b(1) * (Ud(1) - V(1)) - Psub - Q(1); % F(1)

% For 2...M, F(j)=B(j-1)*(Ud(j-1)-V(j-1)) - B(j)*(Ud(j)-V(j)) - Q(j)
for j = 2:M
    F(j)=b(j-1)*(Ud(j-1)-V(j-1)) - b(j)*(Ud(j)-V(j)) - Q(j);
end

    F(M+1) = b(M) * (Ud(M+1) - V(M)) + b(M+1) * V(M) - Q(M+1); % F(M+1)
    
    % Create linear system and solve
Hs = zeros(M+1);
Hs(1,1) = H(1);
Hs(1,2) = b(1);
for j=2:M
    Hs(j, j-1) = b(j-1);
    Hs(j,j) = H(j);
    Hs(j, j+1) = b(j);
end
Hs(M+1, M) = b(M);
Hs(M+1, M+1) = H(M+1);

    %H = [H1, b1; b2, H2]; % These will change like we discussed
    %F = [F1; F2];
    U = inv(Hs) * F;
    Ulips(n) = U(M+1);      % U(2) -> U(M+1)
    
    % glottis
    Q(1) = 2 * (Lg(1) + L(1)) * U(1) - Q(1);      % Q(1)
    P(1) = b(1) * (U(1) - U(2) - Ud(1) + V(1));    % P(1)
    Vc(1) = 2 * C(1) * P(1) - Vc(1);               % Vc(1)
    u3(1) = Yw(1) * (P(1) + Qwl(1) - Qwc(1));         % u3(1)
    Qwl(1) = 2 * Lw(1) * u3(1) - Qwl(1);            % Qwl(1)
    Qwc(1) = 2 * Cw(1) * u3(1) + Qwc(1);
    
for j = 2:M
   Q(j) = 2 * (L(j-1) + L(j))*U(j) - Q(j);
   P(j) = b(j) * (U(j) - U(j+1) - Ud(j) + V(j));
   Vc(j) = 2 * C(j) * P(j) - Vc(j);           
   u3(j) = Yw(j) * (P(j) + Qwl(j) - Qwc(j));     
   Qwl(j) = 2 * Lw(j) * u3(j) - Qwl(j);          
   Qwc(j) = 2 * Cw(j) * u3(j) + Qwc(j);
end
    % for j>1 Q(j) = 2 * (L(j-1) + L(j))U(j) - Q(j)
    %P(1) = b(1) * (U(1) - U(2) - Ud + V(1));  % P(1)
    % for j>1 P(j) = b(J) * (U(j) + U(j+1) - Ud(j) + V(j))
    
    % The rest are the same for j>1 (but with P(j) instead of P1)
    %Vc = 2 * C * P1 - Vc;               % Vc(1)
    %u3 = Yw * (P1 + Qwl - Qwc);         % u3(1)
    %Qwl = 2 * Lw * u3 - Qwl;            % Qwl(1)
    %Qwc = 2 * Cw * u3 + Qwc;            % Qwc(1)
    
    % lips
    P(M+1) = b(M+1) * (U(M+1) + V(M+1));      % P(M+1) ... U(M+1) etc.
    Q(M+1) = 2 * L(M) * U(M+1) - Q(M+1);     % Q(M+1)
    V(M+1) = -2 * Srad * P(M+1) + V(M+1);   % V(M+1)
    
end

signal = diff(Ulips);% (diff(U2) - mean(diff(U2)))/std(diff(U2));

plot(signal(1:500));
soundsc(signal, Fs);