%M: total length of tube
%Used .00001 instead of 0 for areas of the vocal tract to
%avoid dividing by 0 errors (page 211)

function[b] = equations(AreaFunction, Xj, Fs)
    sz = size(AreaFunction);
    N = sz(1);
    P0 = 1.14 * 10^-3;
    c = 3.5 * 10^4;
    mu = 1.86 * 10^-4;
    b = 1.4*10^3;
    m = 1.5;
    k = 3 * 10^5;
    Psub = 8*980.39;
    T = 1/Fs;
    
   [Sj, Lj, Rj, Cj, Lw, Rw, Cw, Yw, bj, Hj, Fj, Ud, Q, V, Vc, Qwl, Qwc, Pj, u3] = deal(zeros(sz));
    
    
    for n = 1:N
        Sj(n) = 2 * AreaFunction(n) * pi^(1/2);
        Lj(n) = P0 * Xj / (2*AreaFunction(n));
        Rj(n) = 4 * pi * mu * Xj / AreaFunction(n);
        Cj(n) = Xj * AreaFunction(n) / ((P0 * c)^2);
        Lw(n) = m / (Xj * (Sj(n))^2);
        Rw(n) = b / (Xj * (Sj(n))^2);
        Cw(n) = Xj * (Sj(n))^2 / k;
        Yw(n) = 1 / (2*Lw(n)/T - Rw(n) + T/2*Cw(n));
        bj(n) = 1 / (2*Cj(n)/T + Yw(n));   
    end
end

