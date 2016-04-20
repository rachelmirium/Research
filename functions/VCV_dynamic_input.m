function [A, x] = VCV_dynamic_input(sound1, sound2, sound3)

v1 = load(sound1);
c = load(sound2);
v2 = load(sound3);

M = 100;

A1 = v1.are_A(1,:);
A2 = c.are_A(1,:);
A3 = v2.are_A(1,:);
X1 = v1.are_x(1,:);
X2 = c.are_x(1,:);
X3 = v2.are_x(1,:);

X1c = cumsum(X1);
X2c = cumsum(X2);
X3c = cumsum(X3);

X1i = linspace(X1c(1), X1c(end), M);
X2i = linspace(X2c(1), X2c(end), M);
X3i = linspace(X3c(1), X3c(end), M);

A1i = interp1(X1c, A1, X1i);
A2i = interp1(X2c, A2, X2i);
A3i = interp1(X3c, A3, X3i);

X1id = diff([0 X1i]);
X2id = diff([0 X2i]);
X3id = diff([0 X3i]);

A = zeros(10000, M);
x = zeros(10000, M);

for j = 1:M
    firstA = linspace(A1i(j), A2i(j), 5000);
    secondA = linspace(A2i(j), A3i(j), 5000);
    newA = [firstA secondA];
    
    firstX = linspace(X1id(j), X2id(j), 5000);
    secondX = linspace(X2id(j), X3id(j), 5000);
    newX = [firstX secondX];
    
    A(:,j) = transpose(newA);
    x(:,j) = transpose(newX);
    
    A = max(A, .1);
end

end

