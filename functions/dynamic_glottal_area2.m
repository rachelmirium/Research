%function [Agp] = dynamic_glottal_area(Ap, F0, Fs)

N = length(Ap);
Agp = zeros(N, 1);
dt = 1.0 / Fs;

T0 = 0.01;%1.0 / F0(100);
t1 = .36 * T0;
t2 = .26 * T0;
a = pi / t1;
b = 1.0 / (1.0 - cos(a * t2));

n1 = 0 ;

for n = 1:N
    
    %t = mod((n * dt), T0);
    t = (n - n1) * dt;
    %pause() 
      
    if t < t1
        Agp(n) = .5 * Ap(n) * (1.0 - cos(a * t));
    elseif t < (t1 + t2)
        Agp(n) = Ap(n) * (1.0 - b + b * cos(a * (t - t1)));
    else
        Agp(n) = .00001;
    end
    
    if t > T0
        
        n1 = n
        T0 = 1.0 / F0(n);
        t1 = .36 * T0;
        t2 = .26 * T0;
        a = pi / t1;
        b = 1.0 / (1.0 - cos(a * t2));
        
    end;
    
    %t = (n - n1) * dt;
    
end



% for n = 1:N
%     T0 = 1.0 / F0(n);
%     t1 = .36 * T0;
%     t2 = .26 * T0;
%     a = pi / t1;
%     b = 1.0 / (1.0 - cos(a * t2));
%     t = mod(((n) * dt), T0);
%     
%     if t < t1
%         Agp(n) = .5 * Ap(n) * (1.0 - cos(a * t));
%     elseif t < (t1 + t2)
%         Agp(n) = Ap(n) * (1.0 - b + b * cos(a * (t - t1)));
%     else
%         Agp(n) = .00001;
%     end
% end

%end

