clear
clc

X = [1 + 2i; 3 - 4i];
Y = [2 - 1i; 1 + 3i];

% Element-wise addition
Z_add = X + Y;

% Element-wise multiplication
Z_mul = X .* Y;

% Element-wise power
Z_pow = X .^ 2;

% Element-wise square root
Z_sqrt = sqrt(Y);

% Display results
disp('Element-wise addition:');
disp(Z_add);

disp('Element-wise multiplication:');
disp(Z_mul);

disp('Element-wise power:');
disp(Z_pow);

disp('Element-wise sqrt:');
disp(Z_sqrt);