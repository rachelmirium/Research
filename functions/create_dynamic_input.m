function [A, x] = create_dynamic_input()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
i = load('../synth_vowels/story_who.mat');
u = load('../synth_vowels/story_heed.mat');

A1 = i.are_A(1,:);
A2 = u.are_A(1,:);
X1 = i.are_x(1,:);
X2 = u.are_x(1,:);

for i = 47:50
    A1(i) = .0001;
    X1(i) = .0001;
end
for i = 43:50
    A2(i) = .0001;
    X2(i) = .0001;
end

A = zeros(10000, 50);
x = zeros(10000, 50);

for j = 1:50
    A(:,j) = transpose(linspace(A1(j), A2(j), 10000));
    x(:,j) = transpose(linspace(X1(j), X2(j), 10000));
end

end

