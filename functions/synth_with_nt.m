function signal = synth_with_nt(A, X, Ag0, Ap, F0, Anc, AN, XN)

N = size(A, 1);
if N ~= size(X, 1) || N ~= length(Ag0) || N ~= length(Ap) || N ~= length(F0)
   error('Array dimensions are mismatched');
end

Nlength = size(AN, 1);

%A = max(A, .1);
M = size(A,2);
maxAp = max(Ap);
Fs = 10000;

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

[Agp] = dynamic_glottal_area(Ap, F0, Fs);
Agp = max(Agp,0.001);
Ag = Ag0 + Agp;

N=length(Ag);

[Ulips, Unose] = deal(zeros(N,1));

Ud = 0; 
Q(1) = 0; Q(2) = 0;
V(1) = 0; V(2) = 0;
Vc(1) = 0; Qwl(1) = 0; Qwc(1) = 0;
[Namp, L, R, C, Lw, Rw, Cw, Yw, b, Ud, H, V, F, U, Q, P, Vc, u3, Qwl, Qwc] = deal(zeros(M+1,1));
[UN, PN, ALN, CN, WRN, WLN, WCN, GWN, BN, HN, RN, UDN, QN, VCN, QWLN, QWCN, VN, FN] = deal(zeros(Nlength+1,1));
U(1)=1;
U(2)=2;
AN1 = AN(1)*XN;

for n=2:N
    [Ac, Xc] = min(A(n,:));
    a1 = kc * rho * (1/(maxAp * maxAp) + 1/(Ac * Ac));
    b1 = 12 * glt_len * glt_len * glt_thck / (maxAp * maxAp * maxAp) + 8 * pi * mu * Xc / (Ac * Ac);
    c1 = -1 * Psub;
    p = [a1 b1 c1];
    r = roots(p);
    Udc = -1*r(1);
   
    Lg = (2/T) * rho * glt_thck / Ag(n);
    Rg = (12 * mu * glt_len^2 * glt_thck) / (Ag(n)^3) + ((1.38/2) * rho * abs(U(1)) / (Ag(n)^2));
    
for j = 2:Nlength
    AX = AN(j)*XN;
    AXinv = 1.0 / AX;
    XDA = XN/AN(j);
    ALN(j) = (1/T) * rho * XDA;
    RN(j) = 4 * pi * mu * AXinv / AN(j);
    CN(j) = (2/T) * AX / (rho * c^2);
    WRN(j) = bval * sqrt(pi) * AXinv;
    WLN(j) = 2 / 5 * m *AXinv / (T * sqrt(pi));
    WCN(j) = T/(2*sqrt(pi)) * k * AXinv;
    GWN(j) = 1.0 / (WRN(j) + WLN(j) + WCN(j));
    BN(j) = 1.0 / (CN(j) + GWN(j));
end

for j = 2:Nlength
    HN(j) = -(ALN(j-1) + ALN(j) + RN(j-1) + RN(j) + BN(j-1) + BN(j));
end

RADSN = T*(1.5 * pi * sqrt(pi * AN(Nlength))) / (8 * rho);
RADGN = 9 * pi^2 * AN(Nlength) / (128 * rho * c);
BN(Nlength+1) = 1.0 / (RADSN + RADGN);
HN(Nlength+1) = -(BN(Nlength) + BN(Nlength+1) + ALN(Nlength));
    
for j = 1:M   
    L(j) = (1/T) * rho * X(n,j) / A(n,j);
    R(j) = 4 * pi * mu  / (X(n,j) * (A(n,j))^2);
    C(j) = (2/T) * X(n,j) * A(n,j) / (rho * c^2);
    Lw(j) = 2 * m / (T * A(n,j) * X(n,j) * sqrt(pi));
    Rw(j) = bval / (A(n,j) * X(n,j) * sqrt(pi));
    Cw(j) = T/(2*sqrt(pi)) * k / (A(n,j) * X(n,j));
    Yw(j) = 1 / (Lw(j) + Cw(j) + Rw(j));
    b(j) = 1 / (C(j) + Yw(j));
    Ud(j) = (A(n,j)*X(n,j) - A(n-1,j)*X(n-1,j)) / T;
    %Namp(j) = Udc * Udc / A(n,j) * 5 * 2 * 10^-8;
end
    Namp(min(Xc+1,M)) = Udc * Udc / A(n,min(Xc+1,M)) * 5 * 2 * 10^-8; %insert noise at segment after constriction
    H(1) = -(Lg + L(1) + Rg + R(1) + b(1));

for j = 2:M
    H(j) = - (L(j-1) + L(j) + R(j-1) + R(j) + b(j-1) + b(j));
end

    %H(K) = -(b(K) + L(K) + R(K));
    %H(K+1) = -(b(K+1) + L(K+1) + R(K+1));
    
    Grad = 9 * pi^2 * A(n,M) / (128 * rho * c);
    Srad = T*(1.5 * pi * sqrt(pi * A(n,M))) / (8 * rho);
    b(M+1) = 1 / (Srad + Grad);
    H(M+1) = -(b(M+1) + b(M) + R(M) + L(M)); 
    
for j = 1:M  
    V(j) = Vc(j) - Yw(j) * (Qwl(j) - Qwc(j));
end

    F(1) = - b(1) * (Ud(1) - V(1)) - Psub - Q(1); 

for j = 2:M
    F(j)=b(j-1)*(Ud(j-1)-V(j-1)) - b(j)*(Ud(j)-V(j)) - Q(j);
end
    %F(K+1) = -B(K+1) * (UD(K+1) - V(K+1)) - Q(K+1);
    %F(K) = B(K) * (UD(K) - V(K)) - Q(K);
    F(M+1) = b(M) * (Ud(M+1) - V(M)) + b(M+1) * V(M) - Q(M+1);

for j = 1:Nlength
    VN(j) = VCN(j) - GWN(j) * (QWLN(j) - QWCN(j));
end
FN(1) = -BN(1) * (UDN(1) - VN(1)) - QN(1);
for j = 2:Nlength
    FN(j) = -BN(j-1)*VN(j-1)+BN(j)*VN(j) - QN(j);
end
FN(Nlength+1) = -BN(Nlength) * VN(Nlength) + BN(Nlength+1) * VN(Nlength+1) - QN(Nlength+1);

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

HNs = zeros(Nlength+1,Nlength+2);
HNs(1,1) = 1;
HNs(1,2) = HN(1);
HNs(1,3) = BN(1);
for j = 2:Nlength
    HNs(j,j) = BN(j-1);
    HNs(j,j+1) = HN(j);
    HNs(j,j+2) = BN(j);
end
HNs(Nlength+1,Nlength+1) = BN(Nlength);
HNs(Nlength+1,Nlength+2) = HN(Nlength+1);

    U = inv(Hs) * F;
    UN = inv(HNs) * FN;
    Ulips(n) = U(M+1);
    Unose(n) = UN(Nlength+1);
    
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

%Q(K) = 2 * L(K) * U(K) - Q(K);
%Q(K+1) = 2*L(K+1) * U(K+1) - Q(K+1);

    % lips
    P(M+1) = b(M+1) * (U(M+1) + V(M+1));
    Q(M+1) = 2 * L(M) * U(M+1) - Q(M+1);
    V(M+1) = -2 * Srad * P(M+1) + V(M+1);
    
QN(1) = 2*ALN(1)*UN(1) - QN(1);
PN(1) = PN(1) - BN(1) * UDN(1);
VCN(1) = 2 * CN(1) * PN(1) - VCN(1);
U3N = GWN(1) * (PN(1)+QWLN(1)-QWCN(1));
QWCN(1) = 2 * WCN(1) * U3N + QWCN(1);
QWLN(1) = 2 * WLN(1) * U3N + QWLN(1);
for j = 2:Nlength
    PN(j) = BN(j)*(UN(j)-UN(j+1) + VN(j));
    QN(j) = 2 * (ALN(j-1) + ALN(j)) * UN(j) - QN(j);
    VCN(j) = 2 * CN(j) * PN(j) - VCN(j);
    U3N = GWN(j) * (PN(j)+QWLN(j)-QWCN(j));
    QWCN(j) = 2 * WCN(j) * U3N + QWCN(j);
    QWLN(j) = 2 * WLN(j) * U3N + QWLN(j);
end

PN(Nlength+1) = BN(Nlength+1) * (UN(Nlength+1) + VN(Nlength+1));
QN(Nlength+1) = 2 * ALN(Nlength) * UN(Nlength+1) - QN(Nlength+1);
VN(Nlength+1) = -2 * RADSN * PN(Nlength+1) + VN(Nlength+1);
    
end
signal = diff(Ulips + Unose);% (diff(U2) - mean(diff(U2)))/std(diff(U2));
plot(signal);
soundsc(signal, Fs);