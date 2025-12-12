%% FINITE ELEMENT ANALYSIS OF THIN PLATE: PLANE STRESS FORMULATION
% Simply supported left and bottom edges with horizontal distributed load
% Developed by: [Your Name/Group Name]
% Date: [Current Date]

clear all; close all; clc;

%% 1. PROBLEM PARAMETERS
E = 200e9;          % Young's modulus [Pa]
nu = 0.3;           % Poisson's ratio [-]
t = 0.02;           % Plate thickness [m]
L = 1.5;            % Plate length in x-direction [m]
H = 0.5;            % Plate height in y-direction [m]
p = 1500;           % Distributed load intensity [N/m]

% Convert to traction
traction_x = p;     % Horizontal rightward load [N/m]

%% 2. MESH GENERATION
nx = 30;            % Elements in x-direction
ny = 10;            % Elements in y-direction
nnodes = (nx+1)*(ny+1);     % Total nodes
nelems = nx*ny;             % Total elements

% Generate nodal coordinates
x = linspace(0, L, nx+1);
y = linspace(0, H, ny+1);
[X, Y] = meshgrid(x, y);
nodes = [X(:), Y(:)];

% Element connectivity (4-node quadrilateral)
elements = zeros(nelems, 4);
for i = 1:ny
    for j = 1:nx
        elem_id = (i-1)*nx + j;
        n1 = (i-1)*(nx+1) + j;
        n2 = n1 + 1;
        n3 = n2 + (nx+1);
        n4 = n3 - 1;
        elements(elem_id, :) = [n1, n2, n3, n4];
    end
end

%% 3. MATERIAL MATRIX (PLANE STRESS)
D = (E/(1-nu^2)) * [1, nu, 0;
                   nu, 1, 0;
                   0, 0, (1-nu)/2];

%% 4. GLOBAL STIFFNESS MATRIX ASSEMBLY
ndof = 2 * nnodes;          % Total degrees of freedom
K = sparse(ndof, ndof);     % Sparse storage
F = zeros(ndof, 1);         % Force vector

% Gauss quadrature (2×2)
gauss_points = [-1/sqrt(3), -1/sqrt(3);
                 1/sqrt(3), -1/sqrt(3);
                 1/sqrt(3),  1/sqrt(3);
                -1/sqrt(3),  1/sqrt(3)];
weights = [1, 1, 1, 1];

% Element loop for stiffness matrix
for e = 1:nelems
    % Element nodes and coordinates
    node_ids = elements(e, :);
    elem_nodes = nodes(node_ids, :);
    
    % Initialize element stiffness matrix
    ke = zeros(8, 8);
    
    % Numerical integration
    for gp = 1:4
        xi = gauss_points(gp, 1);
        eta = gauss_points(gp, 2);
        
        % Shape functions
        N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), ...
                    (1+xi)*(1+eta), (1-xi)*(1+eta)];
        
        % Shape function derivatives
        dN_dxi = 0.25 * [-(1-eta),  (1-eta),  (1+eta), -(1+eta)];
        dN_deta = 0.25 * [-(1-xi), -(1+xi),  (1+xi),  (1-xi)];
        
        % Jacobian matrix
        J = [dN_dxi * elem_nodes(:,1), dN_dxi * elem_nodes(:,2);
             dN_deta * elem_nodes(:,1), dN_deta * elem_nodes(:,2)];
        
        detJ = det(J);
        invJ = inv(J);
        
        % Strain-displacement matrix B
        B = zeros(3, 8);
        for i = 1:4
            dN_dx = invJ(1,1)*dN_dxi(i) + invJ(1,2)*dN_deta(i);
            dN_dy = invJ(2,1)*dN_dxi(i) + invJ(2,2)*dN_deta(i);
            
            B(1, 2*i-1) = dN_dx;
            B(2, 2*i) = dN_dy;
            B(3, 2*i-1) = dN_dy;
            B(3, 2*i) = dN_dx;
        end
        
        % Element stiffness contribution
        ke = ke + B' * D * B * detJ * weights(gp) * t;
    end
    
    % Assembly to global matrix
    dof_indices = zeros(8, 1);
    for i = 1:4
        dof_indices(2*i-1:2*i) = [2*node_ids(i)-1, 2*node_ids(i)];
    end
    
    K(dof_indices, dof_indices) = K(dof_indices, dof_indices) + ke;
end

%% 5. APPLY DISTRIBUTED LOAD (RIGHT EDGE)
right_nodes = find(abs(nodes(:,1) - L) < 1e-6);
[~, idx] = sort(nodes(right_nodes, 2));
right_nodes = right_nodes(idx);

% Calculate equivalent nodal forces
for i = 1:length(right_nodes)
    node_id = right_nodes(i);
    
    % Determine load distribution length
    if i == 1
        dy = nodes(right_nodes(2),2) - nodes(node_id,2);
        dy = dy / 2;
    elseif i == length(right_nodes)
        dy = nodes(node_id,2) - nodes(right_nodes(end-1),2);
        dy = dy / 2;
    else
        dy = (nodes(right_nodes(i+1),2) - nodes(right_nodes(i-1),2)) / 2;
    end
    
    % Horizontal force at node
    F(2*node_id-1) = F(2*node_id-1) + traction_x * dy;
end

%% 6. BOUNDARY CONDITIONS
% Left edge nodes (x=0): constraint u=0
left_nodes = find(abs(nodes(:,1)) < 1e-6);

% Bottom edge nodes (y=0): constraint v=0
bottom_nodes = find(abs(nodes(:,2)) < 1e-6);

constrained_dofs = [];

% Constrain x-displacement on left edge
for i = 1:length(left_nodes)
    constrained_dofs = [constrained_dofs, 2*left_nodes(i)-1];
end

% Constrain y-displacement on bottom edge
for i = 1:length(bottom_nodes)
    constrained_dofs = [constrained_dofs, 2*bottom_nodes(i)];
end

constrained_dofs = unique(constrained_dofs);

%% 7. SOLVE SYSTEM OF EQUATIONS
free_dofs = setdiff(1:ndof, constrained_dofs);
U = zeros(ndof, 1);
U(free_dofs) = K(free_dofs, free_dofs) \ F(free_dofs);

% Extract displacement components
u_x = U(1:2:end);   % x-displacement
u_y = U(2:2:end);   % y-displacement

%% 8. RESULTS SUMMARY
fprintf('============================================\n');
fprintf('FINITE ELEMENT ANALYSIS RESULTS\n');
fprintf('============================================\n');
fprintf('Geometry: L=%.2f m, H=%.2f m, t=%.3f m\n', L, H, t);
fprintf('Material: E=%.0f GPa, ν=%.2f\n', E/1e9, nu);
fprintf('Load: p=%.0f N/m (horizontal rightward)\n', p);
fprintf('Mesh: %d nodes, %d elements (30×10)\n', nnodes, nelems);
fprintf('============================================\n');
fprintf('Displacement Results:\n');
fprintf('Max Ux: %.4e m (%.3f mm)\n', max(u_x), max(u_x)*1000);
fprintf('Min Ux: %.4e m (%.3f mm)\n', min(u_x), min(u_x)*1000);
fprintf('Max Uy: %.4e m (%.3f mm)\n', max(u_y), max(u_y)*1000);
fprintf('Min Uy: %.4e m (%.3f mm)\n', min(u_y), min(u_y)*1000);
fprintf('============================================\n');

% Calculate specific point displacements
ur_node = find(abs(nodes(:,1)-L)<1e-6 & abs(nodes(:,2)-H)<1e-6);
lr_node = find(abs(nodes(:,1)-L)<1e-6 & abs(nodes(:,2))<1e-6);
ul_node = find(abs(nodes(:,1))<1e-6 & abs(nodes(:,2)-H)<1e-6);
center_node = find(abs(nodes(:,1)-L/2)<0.01 & abs(nodes(:,2)-H/2)<0.01);

if ~isempty(ur_node)
    fprintf('Top-right corner: Ux=%.4f mm, Uy=%.4f mm\n', ...
            u_x(ur_node)*1000, u_y(ur_node)*1000);
end
if ~isempty(lr_node)
    fprintf('Bottom-right corner: Ux=%.4f mm, Uy=%.4f mm\n', ...
            u_x(lr_node)*1000, u_y(lr_node)*1000);
end
if ~isempty(ul_node)
    fprintf('Top-left corner: Ux=%.4f mm, Uy=%.4f mm\n', ...
            u_x(ul_node)*1000, u_y(ul_node)*1000);
end
if ~isempty(center_node)
    fprintf('Center point: Ux=%.4f mm, Uy=%.4f mm\n', ...
            u_x(center_node)*1000, u_y(center_node)*1000);
end

%% 9. CREATE 2×3 PLOTS LAYOUT
figure('Position', [50, 50, 1400, 800], 'Name', 'Finite Element Analysis Results');

% Prepare data for plotting
x_grid = reshape(nodes(:,1), [ny+1, nx+1]);
y_grid = reshape(nodes(:,2), [ny+1, nx+1]);
u_x_grid = reshape(u_x, [ny+1, nx+1]) * 1000; % Convert to mm
u_y_grid = reshape(u_y, [ny+1, nx+1]) * 1000; % Convert to mm

% Find nodes along specific paths for line plots
top_edge_nodes = find(abs(nodes(:,2) - H) < 1e-6);
[~, idx] = sort(nodes(top_edge_nodes, 1));
top_edge_nodes = top_edge_nodes(idx);

right_edge_nodes = right_nodes; % Already sorted

center_y_nodes = find(abs(nodes(:,2) - H/2) < 1e-6);
[~, idx] = sort(nodes(center_y_nodes, 1));
center_y_nodes = center_y_nodes(idx);

center_x_nodes = find(abs(nodes(:,1) - L/2) < 1e-6);
[~, idx] = sort(nodes(center_x_nodes, 2));
center_x_nodes = center_x_nodes(idx);

%% PLOT 1: Horizontal Displacement along Top Edge
subplot(2, 3, 1);
plot(nodes(top_edge_nodes, 1), u_x(top_edge_nodes)*1000, 'b-', 'LineWidth', 2);
hold on;
plot(nodes(top_edge_nodes, 1), u_y(top_edge_nodes)*1000, 'r--', 'LineWidth', 2);
grid on;
xlabel('x-coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Displacement [mm]', 'FontSize', 10, 'FontWeight', 'bold');
title('Displacements along Top Edge (y=H)', 'FontSize', 11, 'FontWeight', 'bold');
legend('U_x (horizontal)', 'U_y (vertical)', 'Location', 'best');
xlim([0, L]);
ylim([min([u_x; u_y])*1000*1.1, max([u_x; u_y])*1000*1.1]);

%% PLOT 2: Displacements along Right Edge
subplot(2, 3, 2);
plot(nodes(right_edge_nodes, 2), u_x(right_edge_nodes)*1000, 'b-', 'LineWidth', 2);
hold on;
plot(nodes(right_edge_nodes, 2), u_y(right_edge_nodes)*1000, 'r--', 'LineWidth', 2);
grid on;
xlabel('y-coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Displacement [mm]', 'FontSize', 10, 'FontWeight', 'bold');
title('Displacements along Right Edge (x=L)', 'FontSize', 11, 'FontWeight', 'bold');
legend('U_x (horizontal)', 'U_y (vertical)', 'Location', 'best');
xlim([0, H]);

%% PLOT 3: Horizontal Displacement Contour
subplot(2, 3, 3);
contourf(x_grid, y_grid, u_x_grid, 20, 'LineStyle', 'none');
colorbar;
caxis([min(u_x_grid(:)), max(u_x_grid(:))]);
xlabel('x-coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('y-coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
title('Horizontal Displacement U_x [mm]', 'FontSize', 11, 'FontWeight', 'bold');
axis equal tight;
colormap(jet);

%% PLOT 4: Vertical Displacement Contour
subplot(2, 3, 4);
contourf(x_grid, y_grid, u_y_grid, 20, 'LineStyle', 'none');
colorbar;
caxis([min(u_y_grid(:)), max(u_y_grid(:))]);
xlabel('x-coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('y-coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
title('Vertical Displacement U_y [mm]', 'FontSize', 11, 'FontWeight', 'bold');
axis equal tight;
colormap(jet);

%% PLOT 5: Displacement Vector Field
subplot(2, 3, 5);
% Downsample for clarity
skip = 3;
x_vec = x_grid(1:skip:end, 1:skip:end);
y_vec = y_grid(1:skip:end, 1:skip:end);
u_x_vec = u_x_grid(1:skip:end, 1:skip:end);
u_y_vec = u_y_grid(1:skip:end, 1:skip:end);

scale_factor = 0.05;
quiver(x_vec(:), y_vec(:), u_x_vec(:)*scale_factor, u_y_vec(:)*scale_factor, ...
       0, 'b', 'LineWidth', 1, 'MaxHeadSize', 0.5);
hold on;
plot(nodes(:,1), nodes(:,2), 'k.', 'MarkerSize', 6);
xlabel('x-coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('y-coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
title('Displacement Vector Field', 'FontSize', 11, 'FontWeight', 'bold');
axis equal tight;
grid on;
xlim([0, L]);
ylim([0, H]);

%% PLOT 6: Displacement along Centerlines
subplot(2, 3, 6);
% Centerline at y = H/2
plot(nodes(center_y_nodes, 1), u_x(center_y_nodes)*1000, 'b-', 'LineWidth', 2);
hold on;
plot(nodes(center_y_nodes, 1), u_y(center_y_nodes)*1000, 'r--', 'LineWidth', 2);

% Centerline at x = L/2
plot(nodes(center_x_nodes, 2), u_x(center_x_nodes)*1000, 'b:', 'LineWidth', 2);
plot(nodes(center_x_nodes, 2), u_y(center_x_nodes)*1000, 'r-.', 'LineWidth', 2);

grid on;
xlabel('Coordinate [m]', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Displacement [mm]', 'FontSize', 10, 'FontWeight', 'bold');
title('Displacements along Centerlines', 'FontSize', 11, 'FontWeight', 'bold');
legend('U_x at y=H/2', 'U_y at y=H/2', 'U_x at x=L/2', 'U_y at x=L/2', ...
       'Location', 'best');

% Add overall title
sgtitle('Finite Element Analysis: Thin Plate Deformation (Plane Stress)', ...
        'FontSize', 14, 'FontWeight', 'bold');

%% 10. SAVE RESULTS
% Save figure
saveas(gcf, 'FEM_Results_2x3_Plots.png');
saveas(gcf, 'FEM_Results_2x3_Plots.fig');

% Create and save results table
results_table = table((1:nnodes)', nodes(:,1), nodes(:,2), ...
                     u_x, u_y, u_x*1000, u_y*1000, ...
                     'VariableNames', {'Node', 'X_m', 'Y_m', ...
                     'Ux_m', 'Uy_m', 'Ux_mm', 'Uy_mm'});

writetable(results_table, 'FEM_Displacement_Results.csv');

% Save workspace
save('FEM_Analysis_Workspace.mat', 'nodes', 'u_x', 'u_y', 'elements', ...
     'E', 'nu', 't', 'L', 'H', 'p');

fprintf('\nResults saved:\n');
fprintf('  - FEM_Results_2x3_Plots.png (Figure)\n');
fprintf('  - FEM_Results_2x3_Plots.fig (MATLAB Figure)\n');
fprintf('  - FEM_Displacement_Results.csv (Data Table)\n');
fprintf('  - FEM_Analysis_Workspace.mat (Workspace)\n');
fprintf('============================================\n');

%% 11. ADDITIONAL ANALYSIS
% Calculate strain energy
U_array = [u_x'; u_y'];
U_array = U_array(:);
strain_energy = 0.5 * U_array' * K * U_array;
fprintf('Strain Energy: %.4f J\n', strain_energy);

% Calculate average displacements
avg_u_x = mean(u_x);
avg_u_y = mean(u_y);
fprintf('Average Ux: %.4e m\n', avg_u_x);
fprintf('Average Uy: %.4e m\n', avg_u_y);

% Calculate maximum principal strain (approximate)
max_abs_u_x = max(abs(u_x));
max_abs_u_y = max(abs(u_y));
approx_max_strain_x = max_abs_u_x / L;
approx_max_strain_y = max_abs_u_y / H;
fprintf('Approx. max strain in x-direction: %.4e\n', approx_max_strain_x);
fprintf('Approx. max strain in y-direction: %.4e\n', approx_max_strain_y);