%Computes the theta derivatives of each vertex coordinate for the Jansen linkage
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
% these are assumed to be legal values that are roots of the error funcs!
%leg_params: a struct containing the parameters that describe the linkage
%theta: the current angle of the crank
%OUTPUTS:
%dVdtheta: a column vector containing the theta derivates of each vertex coord
function dVdtheta = compute_velocities(vertex_coords, leg_params, theta)
    %your code here
    num_vertices=leg_params.num_vertices;
    num_linkages=leg_params.num_linkages;
    link_to_vertex_list=leg_params.link_to_vertex_list;
    link_lengths=leg_params.link_lengths;
    crank_length=leg_params.crank_length;
    vertex_pos0=leg_params.vertex_pos0;
    vertex_pos2=leg_params.vertex_pos2;

    vertex_coords_mat=column_to_matrix(vertex_coords);

    error_vec = @(vertex_coords) linkage_error_func(vertex_coords, leg_params, theta);

    M=approximate_jacobian(error_vec,vertex_coords);

    B=[zeros(num_linkages,1);-crank_length*sin(theta);crank_length*cos(theta);0;0];

    dVdtheta=M\B;
end