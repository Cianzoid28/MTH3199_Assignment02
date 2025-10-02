function strandbeest_animation()
    
    
    %initialize leg_params structure
    leg_params = struct();
    %number of vertices in linkage
    leg_params.num_vertices = 7;
    %number of links in linkage
    leg_params.num_linkages = 10;
    %matrix relating links to vertices
    leg_params.link_to_vertex_list = ...
    [ 1, 3;... %link 1 adjacency
    3, 4;... %link 2 adjacency
    2, 3;... %link 3 adjacency
    2, 4;... %link 4 adjacency
    4, 5;... %link 5 adjacency
    2, 6;... %link 6 adjacency
    1, 6;... %link 7 adjacency
    5, 6;... %link 8 adjacency
    5, 7;... %link 9 adjacency
    6, 7 ... %link 10 adjacency
    ];
    %length of crank shaft
    leg_params.crank_length = 15.0;
    %fixed position coords of vertex 0
    leg_params.vertex_pos0 = [0;0];
    %fixed position coords of vertex 2
    leg_params.vertex_pos2 = [-38.0;-7.8];
    
    vertex_coords = [...
     15; 0;... %vertex 1 guess
    -38; -7.8;... %vertex 2 guess
    -50; 50;... %vertex 3 guess
    -100; 0;... %vertex 4 guess
    -100; -50;... %vertex 5 guess
    -50; -50;... %vertex 6 guess
    -50; -100]; %vertex 7 guess
    
    %list of lengths for each link
    %in the leg mechanism
    leg_params.link_lengths = ...
    [ 50.0,... %link 1 length
    55.8,... %link 2 length
    41.5,... %link 3 length
    40.1,... %link 4 length
    39.4,... %link 5 length
    39.3,... %link 6 length
    61.9,... %link 7 length
    36.7,... %link 8 length
    65.7,... %link 9 length
    49.0 ... %link 10 length
    ];
plot
   hold on;
    axis([-120,20,-100,50]);
    
    num_vertices = size(vertex_coords,1);
    num_links    = size(leg_params.link_to_vertex_list,1);
    
    % Fixed reference point
    plot(0,0,'ro','MarkerFaceColor','r','MarkerSize',8);
    
    % Plot moving vertices
    square_plot = plot(zeros(num_vertices,1), zeros(num_vertices,1), 'ko', ...
                       'MarkerFaceColor','k','MarkerSize',6);
    
    % Pre-create line objects for each linkage
    link_plots = gobjects(num_links,1); % preallocate graphics handles
    for i = 1:num_links
        link_plots(i) = plot([0 0],[0 0],'k-','LineWidth',2);
    end

    %dotted line for drive
    origin_link = plot([0 0],[0 0],'r:','LineWidth',1.5);
    
    %path of leg
    path_plot = plot(nan,nan,'b-','LineWidth',1.5); % green trajectory
    x_path = []; 
    y_path = [];
    
    for t = 0:0.002:13
        % compute vertex positions
        linkage_point_column   = compute_coords(vertex_coords, leg_params, t);
        linkage_points = column_to_matrix(linkage_point_column);
        
        x_plot = linkage_points(:,1);
        y_plot = linkage_points(:,2);
        
        % update vertices
        set(square_plot,'XData',x_plot,'YData',y_plot);
        
        % update each linkage using the adjacency list
        for i = 1:num_links
            v1 = leg_params.link_to_vertex_list(i,1);
            v2 = leg_params.link_to_vertex_list(i,2);
            set(link_plots(i), 'XData', [x_plot(v1) x_plot(v2)], 'YData', [y_plot(v1) y_plot(v2)]);
        end
        % Update dotted line from origin (0,0) to vertex 1 (crank arm)
        set(origin_link,'XData',[0, x_plot(1)], 'YData',[0, y_plot(1)]);
        
        %Path tracking
        x_path(end+1) = linkage_points(end, 1);
        y_path(end+1) = linkage_points(end, 2);
        set(path_plot,'XData',x_path,'YData',y_path);

    
    drawnow;
    end
end