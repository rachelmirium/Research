function create_airflow_graph(A, X)
    figure;
    title('Airflow vs supraglottal constriction size');
    xlabel('Cross-sectional area of supraglottal constriction (cm^2)');
    ylabel('Volume velocity (cm^3/s)');
    glottal_areas = [.02 .05 .1 .2 .5];
    hold on;
    points = zeros(10000, 2);
    for i = glottal_areas
        points = create_points(A, X, i);
        Ac = points(:,1);
        Udc = points(:, 2);
        plot(Ac, 1/100*Udc);
    end
    legend('Ag = .5 cm^2', 'Ag = .2 cm^2', 'Ag = .1 cm^2', 'Ag = .05 cm^2', 'Ag = .02 cm^2');
end


function points = create_points(A, X, Ap)
M = size(A,2);
%Ap = 0.2;
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
[Namp, L, R, C, Lw, Rw, Cw, Yw, b, Ud, H, V, F, U, Q, P, Vc, u3, Qwl, Qwc] = deal(zeros(M+1,1));
U(1)=1;
U(2)=2;

[Ac, Udc] = deal(zeros(N, 1));
points = zeros(N, 2);

for n=2:N
    [Ac(n), Xc] = min(A(n,:));
    
    a1 = kc * rho * (1/(Ap * Ap) + 1/(Ac(n) * Ac(n)));
    b1 = 12 * glt_len * glt_len * glt_thck / (Ap * Ap * Ap) + 8 * pi * mu * Xc / (Ac(n) * Ac(n));
    c1 = -1 * Psub;
    p = [a1 b1 c1];
    r = roots(p);
    Udc(n) = -1*r(1);
    
    points(n, 1) = Ac(n);
    points(n, 2) = Udc(n);
   
    Lg = (2/T) * rho * glt_thck / Agp(n);
    Rg = (12 * mu * glt_len^2 * glt_thck) / (Agp(n)^3) + ((1.38/2) * rho * abs(U(1)) / (Agp(n)^2));
    
for j = 1:M   
    L(j) = (1/T) * rho * X(n,j) / A(n,j);
    R(j) = 4 * pi * mu  / (X(n,j) * (A(n,j))^2);
    C(j) = (2/T) * X(n,j) * A(n,j) / (rho * c^2);
    Lw(j) = 2 * m / (T * A(n,j) * X(n,j) * sqrt(pi));
    Rw(j) = bval / (A(n,j) * X(n,j) * sqrt(pi));
    Cw(j) = T/(2*sqrt(pi)) * k / (A(n,j) * X(n,j));
    Yw(j) = 1 / (Lw(j) + Cw(j) + Rw(j));
    b(j) = 1 / (C(j) + Yw(j));
    %Namp(j) = Udc(n) * Udc(n) / A(n,j) * 5 * 2 * 10^-8;
    Ud(j) = (A(n,j)*X(n,j) - A(n-1,j)*X(n-1,j)) / T;
end
    Namp(Xc+1) = Udc(n) * Udc(n) / A(n,Xc+1) * 5 * 2 * 10^-8;
    H(1) = -(Lg + L(1) + Rg + R(1) + b(1));

for j = 2:M
    H(j) = - (L(j-1) + L(j) + R(j-1) + R(j) + b(j-1) + b(j));
end
    
    Grad = 9 * pi^2 * A(n,M) / (128 * rho * c);
    Srad = T*(1.5 * pi * sqrt(pi * A(n,M))) / (8 * rho);
    b(M+1) = 1 / (Srad + Grad);
    
    H(M+1) = -(b(M+1) + b(M) + R(M) + L(M)); 
    
    % refresh force constants
for j = 1:M  
    V(j) = Vc(j) - Yw(j) * (Qwl(j) - Qwc(j));
end

    F(1) = - b(1) * (Ud(1) - V(1)) - Psub - Q(1); 

for j = 2:M
    F(j)=b(j-1)*(Ud(j-1)-V(j-1)) - b(j)*(Ud(j)-V(j)) - Q(j);
end

    F(M+1) = b(M) * (Ud(M+1) - V(M)) + b(M+1) * V(M) - Q(M+1);
    
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


    U = inv(Hs) * F;
    Ulips(n) = U(M+1);
    
    % glottis
    Q(1) = 2 * (Lg(1) + L(1)) * U(1) - Q(1);
    P(1) = b(1) * (U(1) - U(2) - Ud(1) + V(1));
    Vc(1) = 2 * C(1) * P(1) - Vc(1);
    u3(1) = Yw(1) * (P(1) + Qwl(1) - Qwc(1));
    Qwl(1) = 2 * Lw(1) * u3(1) - Qwl(1);
    Qwc(1) = 2 * Cw(1) * u3(1) + Qwc(1);
    
for j = 2:M
   Q(j) = 2 * (L(j-1) + L(j))*U(j) - Q(j);
   P(j) = b(j) * (U(j) - U(j+1) - Ud(j) + V(j)) + rand(1,1) * Namp(j);
   Vc(j) = 2 * C(j) * P(j) - Vc(j);           
   u3(j) = Yw(j) * (P(j) + Qwl(j) - Qwc(j));     
   Qwl(j) = 2 * Lw(j) * u3(j) - Qwl(j);          
   Qwc(j) = 2 * Cw(j) * u3(j) + Qwc(j);
end

    % lips
    P(M+1) = b(M+1) * (U(M+1) + V(M+1));
    Q(M+1) = 2 * L(M) * U(M+1) - Q(M+1);
    V(M+1) = -2 * Srad * P(M+1) + V(M+1);
    
end

%signal = diff(Ulips);% (diff(U2) - mean(diff(U2)))/std(diff(U2));

%plot(signal(1:500));
%soundsc(signal, Fs);

points = sortrows(points);
end