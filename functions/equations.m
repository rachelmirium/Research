%M: total length of tube
%Used .00001 instead of 0 for areas of the vocal tract to
%avoid dividing by 0 errors (page 211)

function[] = equations(M, Ap, F0, Fs, D)
    A = calcA(M, Ap, F0, Fs, D);
    X = calcX(M, Fs, D);
    S = calcS(A);
    
    L = calcL(A, Fs);
    R = calcR(A, Fs);
    C = calcC(A, Fs);
    Lw = calcLw(S, Fs);
    Rw = calcRw(S, Fs);
    Cw = calcCw(S, Fs);
    Yw = calcYw(Lw, Rw, Cw, Fs);
    B = calcB(C, Yw, Fs);
    H = calcH(L, R, B, Fs);
    Ud = calcUd(A, X, Fs);
    
    %Q = calcQ(L, U, Fs); %what is U?
    
    [V, Vc, Qwl, Qwc, P, U3] = calcRest(C, Lw, Cw, Yw, B, U, Ud, Fs); %what is U?
    
    %F = calcF(B, Ud, V, Q);
    
    %old tries
    %V = calcV(Vc, Yw, Qwl, Qwc);
    %Vc = calcVc(C, P, Fs);
    %Qwl = calcQwl(Lw, U3, Fs);
    %Qwc = calcQwc(Cw, U3, Fs);
    %P = calcP(B, U, Ud, V);
    %U3 = calcU3(Yw, P, Qwl, Qwc);
end

function A = calcA(M, Ap, F0, Fs, D)
    A = zeros(M, Fs * D);
    for j = 1:M
        A(j,:) = glottal_area(Ap, F0, Fs ,D);
    end
end

function X = calcX(M, Fs, D)
    X = zeros(M, Fs*D);
    for j = 1:M
        X(j,:) = 1/Fs;
    end
end

function S = calcS(A)
    S = 2 * A * pi^(1/2);
end

function L = calcL(A, Fs)
    sz = size(A);
    L = zeros(sz);
    for j = 1:sz(1)
        for n = 1:sz(2)
            L(j,n) = (1.14 * 10^-3) * (1/Fs) / (2 * A(j,n));
        end
    end
end

function R = calcR(A, Fs)
    sz = size(A);
    R = zeros(sz);
    for j = 1:sz(1)
        for n = 1:sz(2)
            R(j,n) = 4 * pi * (1.86 * 10^-4) * (1/Fs) / A(j,n);
        end
    end
end

function C = calcC(A, Fs)
    C = A * (1/Fs) / ((1.14 * 10^-3) * (3.5 * 10^4)^2);
end

function Lw = calcLw(S, Fs)
    sz = size(S);
    Lw = zeros(sz);
    for j = 1:sz(1)
        for n = 1:sz(2)
            Lw(j,n) = (1.5) / (1/Fs * (S(j,n))^2);
        end
    end 
end

function Rw = calcRw(S, Fs)
    sz = size(S);
    Rw = zeros(sz);
    for j = 1:sz(1)
        for n = 1:sz(2)
            Rw(j,n) = (1.4 * 10^3) / (1/Fs * (S(j,n))^2);
        end
    end 
end

function Cw = calcCw(S, Fs)
    sz = size(S);
    Cw = zeros(sz);
    for j = 1:sz(1)
        for n = 1:sz(2)
            Cw(j,n) = (1/Fs * (S(j,n))^2) / (3 * 10^5);
        end
    end 
end

function B = calcB(C, Yw, Fs)
    sz = size(C);
    B = zeros(sz);
    for j = 1:sz(1)
        for n = 1:sz(2)
            B(j, n) = 1 / (2 * C(j, n) / (1/Fs) + Yw(j, n));
        end
    end
end

function Yw = calcYw(Lw, Rw, Cw, Fs)
    sz = size(Lw);
    Yw = zeros(sz);
    for j = 1:sz(1)
        for n = 1:sz(2)
            Yw(j,n) = 1/(2*Lw(j,n)/(1/Fs)-Rw(j,n)+(1/Fs)*2*Cw(j,n));
        end
    end
end

function H = calcH(L, R, B, Fs)
    sz = size(L);
    H = zeros(sz);
    for j = 2:sz(1)
        for n = 1:sz(2)
            H(j,n) = -2*(L(j-1,n)+L(j,n))/(1/Fs)-R(j-1,n)-R(j,n)-B(j-1,n)-B(j,n);
        end
    end
end

function F = calcF(B, Ud, V, Q)
    sz = size(B);
    F = zeros(sz);
    for n = 1:sz(2)
        F(1,n) = -B(1,n) * (Ud(1,n) - V(1,n)) - 8*980.39 - Q(1,n);
    end
    for j = 2:sz(1)
        for n = 2:sz(2)
            F(j,n) = B(j-1,n)*(Ud(j-1,n)-V(j-1,n-1))-B(j,n)*(Ud(j,n)-V(j,n-1)) - Q(j,n-1);
        end
    end
end

function Ud = calcUd(A, X, Fs)
    sz = size(A);
    Ud = zeros(sz);
    for j = 1:sz(1)
        for n = 2:sz(2)
            Ud(j,n) = (A(j,n)*X(j,n)-A(j,n-1)*X(j,n-1)) / (1/Fs);
        end
    end
end

function Q = calcQ(L, U, Fs)
    sz = size(L);
    Q = zeros(sz);
    for j = 2:sz(1)
        for n = 3:sz(2)
            Q(j,n-1) = (4/(1/Fs))*(L(j-1,n-1)+L(j,n-1))*U(j,n-1)-Q(j,n-2);
        end
    end
end

function [V, Vc, Qwl, Qwc, P, U3] = calcRest(C, Lw, Cw, Yw, B, U, Ud, Fs)
    sz = size(Lw);
    V = zeros(sz);
    Vc = zeros(sz);
    Qwl = zeros(sz);
    Qwc = zeros(sz);
    P = zeros(sz);
    U3 = zeros(sz);
    for j = 2:sz(1)
        for n = 3:sz(2)
            Qwl(j,n-1) = (4/(1/Fs))*Lw(j,n-1)*U3(j,n-1)-Qwl(j,n-2);
            Qwc(j,n-1) = ((1/Fs)/Cw(j,n-1))*U3(j,n-1)+Qwc(j,n-2);
            
            Vc(j,n-1) = (4/(1/Fs))*C(j,n-1)*P(j,n-1)-Vc(j,n-2);
            V(j,n-1) = Vc(j,n-1) - Yw(j,n)*(Qwl(j,n-1)-Qwc(j,n-1));

            P(j,n) = B(j,n)*(U(j,n)-U(j-1,n)-Ud(j,n)+V(j,n-1));
            U3(j,n) = Yw(j,n)*(P(j,n)+Qwl(j,n-1)-Qwc(j,n-1));
        end
    end
end

% function V = calcV(Vc, Yw, Qwl, Qwc)
%     sz = size(Vc);
%     V = zeros(sz);
%     for j = 1:sz(1)
%         for n = 2:sz(2)
%             V(j,n-1) = Vc(j,n-1) - Yw(j,n)*(Qwl(j,n-1)-Qwc(j,n-1));
%         end
%     end
% end
% 
% function Vc = calcVc(C, P, Fs)
%     sz = size(C);
%     Vc = zeros(sz);
%     for j = 1:sz(1)
%         for n = 3:sz(2)
%             Vc(j,n-1) = (4/(1/Fs))*C(j,n-1)*P(j,n-1)-Vc(j,n-2);
%         end
%     end
% end
% 
% function Qwl = calcQwl(Lw, U3, Fs)
%     sz = size(Lw);
%     Qwl = zeros(sz);
%     for j = 1:sz(1)
%         for n = 3:sz(2)
%             Qwl(j,n-1) = (4/(1/Fs))*Lw(j,n-1)*U3(j,n-1)-Qwl(j,n-2);
%         end
%     end
% end
% 
% function Qwc = calcQwc(Cw, U3, Fs)
%     sz = size(Cw);
%     Qwc = zeros(sz);
%     for j = 1:sz(1)
%         for n = 3:sz(2)
%            Qwc(j,n-1) = ((1/Fs)/Cw(j,n-1))*U3(j,n-1)+Qwc(j,n-2);
%         end
%     end
% end
% 
% function P = calcP(B, U, Ud, V)
%     sz = size(B);
%     P = zeros(sz);
%     for j = 2:sz(1)
%         for n = 2:sz(2)
%             P(j,n) = B(j,n)*(U(j,n)-U(j-1,n)-Ud(j,n)+V(j,n-1));
%         end
%     end
% end
% 
% function U3 = calcU3(Yw, P, Qwl, Qwc)
%     sz = size(Yw);
%     U3 = zeros(sz);
%     for j = 1:sz(1)
%         for n = 2:sz(2)
%             U3(j,n) = Yw(j,n)*(P(j,n)+Qwl(j,n-1)-Qwc(j,n-1));
%         end
%     end
% end