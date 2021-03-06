function Y = calc_U2(AreaFunction, Fs)

P0 = 1.14 * 10^-3;
c = 3.5 * 10^4;
Psub = 8 * 980.39;
A1 = 17;
X1 = 1;
[L, R, Ud, V, Q, b] = equations(AreaFunction, X1, Fs);

b1 = b(1);
H1 = -(L(1) + L(1) + R(1) + R(1) + b(1));
Srad = 9 * pi * A1 / (128 * P0 * c);
Grad = (3 * pi * (pi * A1)^(1/2)) / (8 * P0);
b2 = 1 / (Srad + Grad);
H2 = -(b(2) + b(1) + L(1) + R(1));
F1 = -b(1) * (Ud(1) - V(1)) - Psub - Q(1);
F2 = b(1) * (Ud(1) - V(1)) + b(2) * V(2) - Q(2);
A = [H1, b1; b2, H2];
Y = [F1; F2];
X = inv(A) * Y;
U2 = X(2,:);
sound(diff(U2), Fs);
end

function[L, R, Ud, V, Q, b] = equations(AreaFunction, Xj, Fs)
    sz = size(AreaFunction);
    N = sz(1);
    P0 = 1.14 * 10^-3;
    c = 3.5 * 10^4;
    mu = 1.86 * 10^-4;
    bval = 1.4*10^3;
    m = 1.5;
    k = 3 * 10^5;
    Psub = 8*980.39;
    T = 1/Fs;
    
    [S, L, R, C, Lw, Rw, Cw, Yw, b, H, F, Ud, Q, V, Vc, Qwl, Qwc, P, u3] = deal(ones(sz));  
    
    for n = 1:N
        S(n) = 2 * AreaFunction(n) * pi^(1/2);
        L(n) = P0 * Xj / (2*AreaFunction(n));
        R(n) = 4 * pi * mu * Xj / AreaFunction(n);
        C(n) = Xj * AreaFunction(n) / ((P0 * c)^2);
        Lw(n) = m / (Xj * (S(n))^2);
        Rw(n) = bval / (Xj * (S(n))^2);
        Cw(n) = Xj * (S(n))^2 / k;
        Yw(n) = 1 / (2*Lw(n)/T - Rw(n) + T/2*Cw(n));
        b(n) = 1 / (2*C(n)/T + Yw(n)); 
        H(n) = -2*(L(n) + L(n))/T - R(n) - R(n) - b(n) - b(n);
        Ud(n) = (AreaFunction(n) * Xj) / T;
        
        if(n == 2)
            %Q(n-1) = (4/T)*(L(n-1)+L(n-1)) * U(n-1);
            Vc(n-1) = (4/T)*C(n-1)*P(n-1);
            Qwl(n-1) = (4/T)*Lw(n-1)*u3(n-1);
            Qwc(n-1) = (T/Cw(n-1))*u3(n-1);
        end
        
        if(n > 2)
            %Q(n-1) = (4/T)*(L(n-1)+L(n-1)) * U(n-1)-Q(n-2);
            Vc(n-1) = (4/T)*C(n-1)*P(n-1)-Vc(n-2);
            Qwl(n-1) = (4/T)*Lw(n-1)*u3(n-1) - Qwl(n-2);
            Qwc(n-1) = (T/Cw(n-1))*u3(n-1) + Qwc(n-2);
            V(n-1) = Vc(n-1) - Yw(n)*(Qwl(n-1) - Qwc(n-1));
            F(n) = b(n) * (Ud(n) - V(n-1)) - b(n)*(Ud(n)-V(n-1)) - Q(n-1);
            %P(n) = b(n)*(U(n) - U(n) - Ud(n) + V(n-1));
            u3(n) = Yw(n) * (P(n) + Qwl(n-1) - Qwc(n-1));
        end
    end
end