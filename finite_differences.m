%Finite differences method for calculating tip velocity
%vertext_cords: a column vector containing the (x,y) coordinates of every vertex
%leg_params: a struct containing the parameters that describe the linkage
%theta: the current angle of the crank
%OUTPUTS:
%dVdtheta: a column vector containing the theta derivates of each vertex coord
function dvdtheta = finite_differences(vertex_coords, leg_params, theta)
    vertex_coords_root= @(theta) compute_coords(vertex_coords, leg_params, theta);
    
    dvdtheta=approximate_derivative(vertex_coords_root, theta, 1);
    
end