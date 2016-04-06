function [A, x] = create_dynamic_input()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
u = load('../synth_vowels/story_who.mat');
i = load('../synth_vowels/story_heed.mat');

A1 = i.are_A(1,:);
A2 = u.are_A(1,:);
X1 = i.are_x(1,:);
X2 = u.are_x(1,:);

%X1c = [0 cumsum(X1)];
X1c = cumsum(X1);
%X2c = [0 cumsum(X2)];
X2c = cumsum(X2);

%X1s=linspace(.3968,X2(end),46);
X1i = linspace(X1c(1), X1c(end), 100);
X2i = linspace(X2c(1), X2c(end), 100);

%A1 = interp1(X1s, A2, X2);
% size(A1)
% size(X1c)
% size(X1i)
A1i = interp1(X1c, A1, X1i);
A2i = interp1(X2c, A2, X2i);

X1id = diff([0 X1i]);
X2id = diff([0 X2i]);

A = zeros(10000, 100);
x = zeros(10000, 100);

for j = 1:100
    A(:,j) = transpose(linspace(A1i(j), A2i(j), 10000));
    x(:,j) = transpose(linspace(X1id(j), X2id(j), 10000));
end

end

