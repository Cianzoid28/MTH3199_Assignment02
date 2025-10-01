function length_errors = link_length_error_func(vertex_coords, leg_params)

num_links = size(leg_params.link_to_vertex_list, 1);
length_errors = zeros(num_links, 1);

    for i = 1:num_links
    vertexA = vertex_coords(leg_params.link_to_vertex_list(i, 1), :);
    vertexB = vertex_coords(leg_params.link_to_vertex_list(i, 2), :);
    d = leg_params.link_lengths(i);
    length_errors(i) = norm(vertexB - vertexA)^2 - d^2;

    end
end