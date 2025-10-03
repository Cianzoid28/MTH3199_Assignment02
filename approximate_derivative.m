%example of how to implement finite difference approximation
%for the first and second derivative of a function
%INPUTS:
%fun: the mathetmatical function we want to differentiate
%x: the input value of fun that we want to compute the derivative at
%OUTPUTS:
%dfdx: approximation of fun'(x)
function dfdx = approximate_derivative(fun,x,i)

    %set the step size to be tiny
    delta_x = 1e-6;
    ei=zeros(length(x),1);
    ei(i)=1;

    %compute the function at different points near x
    f_left = fun(x-delta_x*ei);
    f_right = fun(x+delta_x*ei);

    %approximate the first derivative
    dfdx = (f_right-f_left)/(2*delta_x);
end