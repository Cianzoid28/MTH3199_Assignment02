%Computes the vertex coordinates that describe a legal linkage configuration
%INPUTS:
%vertex_coords_guess: a column vector containing the (x,y) coordinates of every vertex
% these coords are just a GUESS! It's used to seed Newton's method
%leg_params: a struct containing the parameters that describe the linkage
%theta: the desired angle of the crank
%OUTPUTS:
%vertex_coords_root: a column vector containing the (x,y) coordinates of every vertex
% these coords satisfy all the kinematic constraints!
function vertex_coords_root = compute_coords(vertex_coords_guess, leg_params, theta)
    solver_params = struct();
    solver_params.ftol = 1e-10;
    solver_params.dxmax = 1e8;
    solver_params.dxmin = 1e-10;
    solver_params.max_iter = 10;
    solver_params.approx = 1;
    

    error_vec = @(vertex_coords) linkage_error_func(vertex_coords, leg_params, theta);
    vertex_coords_root = multivarate_newton(error_vec, vertex_coords_guess, solver_params);

end