function [A, x] = create_dynamic_input()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
u = load('../synth_vowels/story_who.mat');
i = load('../synth_vowels/story_heed.mat');

A1 = i.are_A(1,:);
A2 = u.are_A(1,:);
X1 = i.are_x(1,:);
X2 = u.are_x(1,:);

X1s=linspace(.3968,X2(end),46);
A1 = interp1(X1s, A2, X2);

A = zeros(10000, 46);
x = zeros(10000, 46);

for j = 1:46
    A(:,j) = transpose(linspace(A1(j), A2(j), 10000));
    x(:,j) = transpose(linspace(X1s(j), X2(j), 10000));
end

end

