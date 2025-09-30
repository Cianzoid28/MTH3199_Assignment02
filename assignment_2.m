function assignment_2()
       %the following is a function made to test our analitical or
       %numerical jacobian during the root finding process.


       X0 = randn(3,1); % Creates 3 random numbers from a standard normal distribution
       [~, J_analytical] = test_function01(X0); % simply creates analytical jacobian from test_function01
       J_numerical = approximate_jacobian(@test_function01, X0); % makes numerical jacobian from approximate_jacobian function

       solver_params= struct(); % makes structure for parameter values
       solver_params.dxmin = 1e-10; % min step size
       solver_params.ftol = 1e-10 ; %tolerence on function value
       solver_params.dxmax  = 1e8; % max step size
       solver_params.max_iter = 200 ; % max # of iterations
       solver_params.approx = 1; % 1 = use approximate jacobian, 0 = use numerical jacobian

       x_root = mulivarate_newton(@test_function01,X0,solver_params); % finds root based on multivariate_newton function using
       % our given function, our initial guesses, and our parameters

       disp(x_root) % displays the found root for out function

       f_root = test_function01(x_root); % plugs in the root into our function to verify if correct
       disp(f_root) % displays the value of the function at the root

       % disp(J_analytical)
       % disp(J_numerical)

end




function [fval, J] = test_function01(x)

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

%% 

%Implementation of finite difference approximation
%for Jacobian of multidimensional function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%J: approximation of Jacobian of fun at x

function J = approximate_jacobian(fun,x)
    f0 = fun(x); % Vector of function with x inputs
    J = zeros(length(f0), length(x)); %makes a matrix of zeros with the length of the vector f0 by the length x

    e_n = zeros(length(x),1); % makes a column vector size x by 1 (so 3 by 1)
    
    delta_x = 1e-6; %Step size

    for n = 1:length(x) % Loop for each value of vector x
        e_n(n) = 1; % Makes standard basis vector

        % The following code is derived from the numerical diffrentiation
        % function for the approximate derivative where: 

        f_left = fun(x - e_n*delta_x); % is the left part of the numerator of the function
        f_right = fun(x + e_n*delta_x); % is the right part of the numeration of the function
        J(:,n) = (f_right - f_left)/(2*delta_x); % the function itself for deriving the approximate derivative

        e_n(n) = 0; % Resetting to 0 so that the vector does not become all ones
       
    end
end

%% 

% This function calculates the root of a function using the
% multidimensional version of newtons method. It does so given a function,
% an input value x from which we want to begin the process from, and a set
% of specified parameters, solver_params


function X_root  = mulivarate_newton(fun,x,solver_params)

    dxmin = solver_params.dxmin; % min step size from specified params from solver_params
    ftol = solver_params.ftol; %tolerence on function value from specified params from solver_params
    dxmax = solver_params.dxmax; % max step size from specified params from solver_params
    max_iter = solver_params.max_iter; % max # of iterations from specified params from solver_params
    approx = solver_params.approx; % 1 = use approximate jacobian, 0 = use numerical jacobian from specified params from solver_params
 
    if approx % if approx is 1, will use approximate jacobian, 0 will use numerical jacobian
        fval = fun(x); % function value at guess x 
        J = approximate_jacobian(fun,x); %uses approx jacobian
    else 
        [fval, J] = fun(x); % uses numerical jacobian
    end 
    delta_x = -J\fval ; % finds step using jacobian

    count = 0; % starts keeping count of iterations

    while count < max_iter && norm(delta_x) > dxmin && norm(fval) > ftol && norm(delta_x) < dxmax % makes sure the max amount of iterations is not passed
% and the guess is greater than the min guess size and value of our guess is greater then out tolerance, and the guess is less than the maximum guess... 
        count = count + 1; % adds to iteration count

        if approx % if approx is 1, will use approximate jacobian, 0 will use numerical jacobian
            fval = fun(x); % function value at guess x 
            J = approximate_jacobian(fun, x); %uses approx jacobian
        else
            [fval, J] = fun(x); % uses numerical jacobian
        end
        delta_x = -J\fval; % finds step using jacobian
        x = x + delta_x; % updates root
    end

    X_root = x; %calls root X_root

end




%% 

function test_numerical_jacobian()
%number of tests to perform
    num_tests = 100;
%iterate num_tests times
        for n = 1:num_tests
%generate a randomized input and output dimension
        input_dim = randi([1,15]);
        output_dim = randi([1,15]);
%generate a input_dim x input_dim x output_dim matrix stack A
        A = randn(input_dim,input_dim,output_dim);
%generate a matrix, B of dimension output_dim x input_dim
        B = randn(output_dim,input_dim);
%generate a column vector, C of height output_dim
        C = randn(output_dim,1);
%create a new test function
%this is essentially a random second-order (quadratic) function
%with input dimension input_dim and output dimension output_dim
        test_fun = @(X) jacobian_test_function(X,A,B,C);
        X_guess = randn(input_dim,1);
        %evaluate numerical Jacobian of test_fun
        %use whatever your function name was here!
        J_numerical = approximate_jacobian(test_fun,X_guess);
        %compute the analytical jacobian of jacobian_test_function
        J_analytical = B;
        for n = 1:output_dim
        J_analytical(n,:)=J_analytical(n,:)+X_guess'*A(:,:,n);
        J_analytical(n,:)=J_analytical(n,:)+X_guess'*A(:,:,n)';
        end
        %compare with Jacobian of A
        largest_error = max(max(abs(J_numerical-J_analytical)));
        %if J is not close to A, print fail.
        if largest_error>1e-7
        disp('fail!');
        end
    end
end


%computes a quadratic function on input X
function f_val = jacobian_test_function(X,A,B,C)
    output_length = length(C);
    f_val = B*X+C;
    for n = 1:output_length
    f_val(n)=f_val(n)+(X'*A(:,:,n)*X);
    end
end





    




