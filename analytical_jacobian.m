function [fval, J] = analytical_jacobian(x)

% The following function stores the function we are working with along with
% the analytical jacobian given a vector x

    x1 = x(1); x2 = x(2); x3= x(3); % established x1 x2 and x3

    f1 = x1^2 + x2^2 - 6 - x3^5; % Establishes the functions for the matrix
    f2 = x1*(x3) + x2 - 12; % ^
    f3 = sin(x1 +x2 + x3); % ^

    fval = [f1;f2;f3]; % makes the actual matrix

    J = [[2*x1, 2*x2, + 5*x3^4];...
        [x3, 1, x1];...
        [cos(x1 + x2 + x3),cos(x1 + x2 + x3),cos(x1 + x2 + x3)]]; % the following is the analytical jacobian through partial diffrentiation
end 