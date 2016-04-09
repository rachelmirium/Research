function [A, x] = create_dynamic_input(sound1, sound2)

s1 = load(sound1);
s2 = load(sound2);

M = 100;

A1 = s1.are_A(1,:);
A2 = s2.are_A(1,:);
X1 = s1.are_x(1,:);
X2 = s2.are_x(1,:);

X1c = cumsum(X1);
X2c = cumsum(X2);

X1i = linspace(X1c(1), X1c(end), M);
X2i = linspace(X2c(1), X2c(end), M);

A1i = interp1(X1c, A1, X1i);
A2i = interp1(X2c, A2, X2i);

X1id = diff([0 X1i]);
X2id = diff([0 X2i]);

A = zeros(10000, M);
x = zeros(10000, M);

for j = 1:M
    A(:,j) = transpose(linspace(A1i(j), A2i(j), 10000));
    x(:,j) = transpose(linspace(X1id(j), X2id(j), 10000));
end

end

