function [A, X, Ag0, Ap, F0] = ASA_inputs()

%new area function 
M = 100;

V = load('../synth_vowels/story_had.mat');
S = load('../synth_vowels/story_cut.mat');

A1 = V.are_A(1,:);
A2 = S.are_A(1,:);
X1 = V.are_x(1,:);
X2 = S.are_x(1,:);

X1c = cumsum(X1);
X2c = cumsum(X2);

X1i = linspace(X1c(1), X1c(end), M);
X2i = linspace(X2c(1), X2c(end), M);

A1i = interp1(X1c, A1, X1i);
A2i = interp1(X2c, A2, X2i);

X1id = diff([0 X1i]);
X2id = diff([0 X2i]);

[A, X] = deal(zeros(4000, M));

for j = 1:M
    x = [1 501 1001 2001 2501 4000];
    Av = [A1i(j) A1i(j) A2i(j) A2i(j) A1i(j) A2i(j)];
    Xv = [X1id(j) X1id(j) X2id(j) X2id(j) X1id(j) X1id(j)];
    
    A(1:500, j) = A1i(j);
    A(501:1000, j) = interp1(x, Av, 501:1000, 'pchip');
    A(1001:2000, j) = A2i(j);
    A(2001:2500, j) = interp1(x, Av, 2001:2500, 'pchip');
    A(2501:4000, j) = A1i(j);
    
    X(1:500, j) = X1id(j);
    X(501:1000, j) = interp1(x, Xv, 501:1000, 'pchip');
    X(1001:2000, j) = X2id(j);
    X(2001:2500, j) = interp1(x, Xv, 2001:2500, 'pchip');
    X(2501:4000, j) = X1id(j);
    
end

%dynamic Ag0, Ap, F0

[Ag0, Ap, F0] = deal(zeros(4000, 1));

x = [1 501 1001 2001 2501 3501 4001];
v = [0 0 .4 .4 0 0 .1];

Ag0(1:500) = 0;
Ag0(501:1000) = interp1(x, v, 501:1000, 'pchip');
Ag0(1001:2000) = .4;
Ag0(2001:2500) = interp1(x, v, 2001:2500, 'pchip');
Ag0(2501:3500) = 0;
Ag0(3501:4000) = interp1(x, v, 3501:4000, 'pchip');

Ap(1:200) = linspace(.1, .2, 200);
Ap(201:500) = .2;
Ap(501:700) = linspace(.2, .1, 200);
Ap(701:2300) = .1;
Ap(2301:2500) = linspace(.1, .2, 200);
Ap(2501:3500) = linspace(.2, .1, 1000);
Ap(3501:4000) = linspace(.1, 0, 500);

F0(1:700) = linspace(100, 96, 700);
F0(701:2300) = 130;
F0(2301:3500) = 100;%linspace(130, 90, 1200);
F0(3501:4000) = linspace(90, 70, 500);

end

