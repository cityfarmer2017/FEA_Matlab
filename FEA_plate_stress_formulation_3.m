%% Finite Element Analysis of Thin Plate: Plane Stress Formulation
% Developed by: [Name(s)]
% Date: [Current Date]
% Problem: Simply supported left and bottom edges with horizontal (direction: left to right) distributed load

clear all; close all; clc;

%% 1. Problem Parameters Definition
E = 200e9;          % Young's modulus (Pa)
nu = 0.3;           % Poisson's ratio
t = 0.02;           % Plate thickness (m)
L = 1.5;            % Plate length in x-direction (m)
H = 0.5;            % Plate height in y-direction (m)
p = 1500;           % Distributed load intensity (N/m) - horizontal from left to right

% Convert to traction (N/m)
traction_x = p;     % Positive indicates rightward direction

%% 2. Mesh Generation
nx = 30;            % Number of elements in x-direction
ny = 10;            % Number of elements in y-direction
nnodes = (nx+1)*(ny+1);     % Total number of nodes
nelems = nx*ny;             % Total number of elements

% Generate element connectivity (4-node quadrilateral elements)
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

% Generate nodal coordinates
[X, Y] = meshgrid(linspace(0, L, nx+1), linspace(0, H, ny+1));
nodes = [X(:), Y(:)];

%% 3. Material Matrix (Plane Stress)
D = (E/(1-nu^2)) * [1, nu, 0;
                    nu, 1, 0;
                    0, 0, (1-nu)/2];

%% 4. Global Stiffness Matrix and Force Vector Assembly
ndof = 2 * nnodes;          % Total degrees of freedom
K = sparse(ndof, ndof);     % Sparse storage for efficiency
F = zeros(ndof, 1);

% Gauss quadrature points (2×2)
gauss_points = [-1/sqrt(3), -1/sqrt(3);
                 1/sqrt(3), -1/sqrt(3);
                 1/sqrt(3),  1/sqrt(3);
                -1/sqrt(3),  1/sqrt(3)];
weights = [1, 1, 1, 1];

% Element loop for stiffness matrix assembly
for e = 1:nelems
    % Get node coordinates for current element
    node_ids = elements(e, :);
    elem_nodes = nodes(node_ids, :);

    % Initialize element stiffness matrix
    ke = zeros(8, 8);

    % Numerical integration using Gauss quadrature
    for gp = 1:4
        xi = gauss_points(gp, 1);
        eta = gauss_points(gp, 2);

        % Shape functions and derivatives
        N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];

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

        % Element stiffness matrix contribution
        ke = ke + B' * D * B * detJ * weights(gp) * t;
    end

    % Assemble into global stiffness matrix
    dof_indices = zeros(8, 1);
    for i = 1:4
        dof_indices(2*i-1:2*i) = [2*node_ids(i)-1, 2*node_ids(i)];
    end
    
    K(dof_indices, dof_indices) = K(dof_indices, dof_indices) + ke;
end

%% 5. Apply Distributed Load (Right Edge - Horizontal Direction)
% Find nodes on right edge (x = L)
right_nodes = find(abs(nodes(:,1) - L) < 1e-6);

% Sort by y-coordinate
[~, idx] = sort(nodes(right_nodes, 2));
right_nodes = right_nodes(idx);

% Calculate equivalent nodal forces
for i = 1:length(right_nodes)
    node_id = right_nodes(i);

    % Determine load distribution length for this node
    if i == 1
        % First node: consider only upper half
        dy = nodes(right_nodes(2),2) - nodes(node_id,2);
        dy = dy / 2;
    elseif i == length(right_nodes)
        % Last node: consider only lower half
        dy = nodes(node_id,2) - nodes(right_nodes(end-1),2);
        dy = dy / 2;
    else
        % Middle nodes: consider both halves
        dy = (nodes(right_nodes(i+1),2) - nodes(right_nodes(i-1),2)) / 2;
    end

    % Equivalent nodal force (horizontal rightward)
    F(2*node_id-1) = F(2*node_id-1) + traction_x * dy;
end

%% 6. Apply Boundary Conditions (Simply Supported sides: Left and Bottom Edges)
% Find nodes on left edge (x = 0)
left_nodes = find(abs(nodes(:,1)) < 1e-6);

% Find nodes on bottom edge (y = 0)
bottom_nodes = find(abs(nodes(:,2)) < 1e-6);

% Constrain degrees of freedom:
% Left edge: constrain horizontal displacement (u=0)
% Bottom edge: constrain vertical displacement (v=0)
constrained_dofs = [];

% Left edge nodes: constrain x-direction displacement (odd-numbered DOFs)
for i = 1:length(left_nodes)
    constrained_dofs = [constrained_dofs, 2*left_nodes(i)-1];
end

% Bottom edge nodes: constrain y-direction displacement (even-numbered DOFs)
for i = 1:length(bottom_nodes)
    constrained_dofs = [constrained_dofs, 2*bottom_nodes(i)];
end

% Remove duplicates (bottom-left corner node belongs to both edges)
constrained_dofs = unique(constrained_dofs);

%% 7. Solve System of Equations
free_dofs = setdiff(1:ndof, constrained_dofs);

% Solve using reduced system approach
U = zeros(ndof, 1);
U(free_dofs) = K(free_dofs, free_dofs) \ F(free_dofs);

% Extract displacement components
u_x = U(1:2:end);   % x-direction displacement
u_y = U(2:2:end);   % y-direction displacement

%% 8. Results Analysis and Visualization
fprintf('=== FINITE ELEMENT ANALYSIS RESULTS ===\n');
fprintf('Total nodes: %d\n', nnodes);
fprintf('Total elements: %d\n', nelems);
fprintf('Maximum horizontal displacement: %.6e m\n', max(u_x));
fprintf('Minimum horizontal displacement: %.6e m\n', min(u_x));
fprintf('Maximum vertical displacement: %.6e m\n', max(u_y));
fprintf('Minimum vertical displacement: %.6e m\n', min(u_y));
fprintf('Bottom-right corner horizontal displacement: %.6e m\n', u_x(right_nodes(end)));
fprintf('Top-right corner horizontal displacement: %.6e m\n', u_x(right_nodes(1)));

% Create results table
results_table = table((1:nnodes)', nodes(:,1), nodes(:,2), u_x, u_y, 'VariableNames', {'Node', 'X', 'Y', 'Ux', 'Uy'});

% Display first 20 nodes results
fprintf('\nDisplacement results for first 20 nodes:\n');
disp(results_table(1:min(20, nnodes), :));

% Calculate displacements along specific paths
% Top edge (y=H)
top_edge_nodes = find(abs(nodes(:,2) - H) < 1e-6);
[~, idx] = sort(nodes(top_edge_nodes, 1));
top_edge_nodes = top_edge_nodes(idx);

% Right edge (x=L)
[~, idx] = sort(nodes(right_nodes, 2));
right_edge_nodes = right_nodes(idx);

%% 9. Visualization of Results
figure('Position', [50, 50, 1200, 900]);

% Subplot 1: Displacements along top edge
subplot(3, 2, 1);
plot(nodes(top_edge_nodes, 1), u_x(top_edge_nodes)*1e3, 'b-o', 'LineWidth', 2, 'DisplayName', 'U_x');
hold on;
plot(nodes(top_edge_nodes, 1), u_y(top_edge_nodes)*1e3, 'r-s', 'LineWidth', 2, 'DisplayName', 'U_y');
xlabel('x-coordinate (m)');
ylabel('Displacement (mm)');
title('Displacements along Top Edge (y=H)');
legend('Location', 'best');
grid on;
xlim([0, L]);

% Subplot 2: Displacements along right edge
subplot(3, 2, 2);
plot(nodes(right_edge_nodes, 2), u_x(right_edge_nodes)*1e3, 'b-o', 'LineWidth', 2, 'DisplayName', 'U_x');
hold on;
plot(nodes(right_edge_nodes, 2), u_y(right_edge_nodes)*1e3, 'r-s', 'LineWidth', 2, 'DisplayName', 'U_y');
xlabel('y-coordinate (m)');
ylabel('Displacement (mm)');
title('Displacements along Right Edge (x=L)');
legend('Location', 'best');
grid on;
xlim([0, H]);

% Subplot 3: Horizontal displacement contour
subplot(3, 2, 3);
% Reshape to grid format for contour plotting
x_grid = reshape(nodes(:,1), [ny+1, nx+1]);
y_grid = reshape(nodes(:,2), [ny+1, nx+1]);
u_x_grid = reshape(u_x, [ny+1, nx+1]);
contourf(x_grid, y_grid, u_x_grid*1e3, 20, 'LineStyle', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Horizontal Displacement U_x (mm)');
axis equal tight;

% Subplot 4: Vertical displacement contour
subplot(3, 2, 4);
u_y_grid = reshape(u_y, [ny+1, nx+1]);
contourf(x_grid, y_grid, u_y_grid*1e3, 20, 'LineStyle', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Vertical Displacement U_y (mm)');
axis equal tight;

% Subplot 5: Displacement vector plot
subplot(3, 2, 5);
scale_factor = 0.1 / max(sqrt(u_x.^2 + u_y.^2)); % Scale for visualization
quiver(nodes(:,1), nodes(:,2), u_x*scale_factor, u_y*scale_factor, 0, 'b');
hold on;
plot(nodes(:,1), nodes(:,2), 'k.', 'MarkerSize', 8);
xlabel('x (m)');
ylabel('y (m)');
title('Displacement Vector Plot (Scaled)');
axis equal tight;
grid on;

% Subplot 6: Deformed vs. undeformed mesh comparison
subplot(3, 2, 6);
% Undeformed mesh
for e = 1:nelems
    node_ids = elements(e, :);
    patch(nodes(node_ids, 1), nodes(node_ids, 2), 'w', 'EdgeColor', [0.7 0.7 0.7], 'FaceColor', 'none');
end
hold on;
% Deformed mesh (scaled for visualization)
deformation_scale = 100; % Magnification factor
deformed_nodes = nodes + [u_x, u_y] * deformation_scale;
for e = 1:min(nelems, 200) % Plot subset to avoid overcrowding
    node_ids = elements(e, :);
    patch(deformed_nodes(node_ids, 1), deformed_nodes(node_ids, 2), 'w', 'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 1.5);
end
xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Mesh Comparison (Deformation Magnified %d times)', deformation_scale));
axis equal tight;
grid on;

% Save results to files
writetable(results_table, 'displacement_results_horizontal_load.csv');
save('FEM_results_horizontal_load.mat', 'nodes', 'u_x', 'u_y', 'elements', 'right_nodes', 'top_edge_nodes');

fprintf('\nResults saved to files:\n');
fprintf('  - displacement_results_horizontal_load.csv (Displacement data table)\n');
fprintf('  - FEM_results_horizontal_load.mat (MATLAB data file)\n');

%% 10. Additional Statistical Analysis
fprintf('\n=== DISPLACEMENT STATISTICS ===\n');
fprintf('Horizontal displacement range: %.4f to %.4f mm\n', min(u_x)*1e3, max(u_x)*1e3);
fprintf('Vertical displacement range: %.4f to %.4f mm\n', min(u_y)*1e3, max(u_y)*1e3);
fprintf('Average horizontal displacement: %.4f mm\n', mean(u_x)*1e3);
fprintf('Average vertical displacement: %.4f mm\n', mean(u_y)*1e3);

% Calculate displacements at specific critical points
% Top-right corner
ur_node = find(abs(nodes(:,1)-L)<1e-6 & abs(nodes(:,2)-H)<1e-6);
if ~isempty(ur_node)
    fprintf('Top-right corner: Ux=%.4f mm, Uy=%.4f mm\n', u_x(ur_node)*1e3, u_y(ur_node)*1e3);
end

% Bottom-right corner
lr_node = find(abs(nodes(:,1)-L)<1e-6 & abs(nodes(:,2))<1e-6);
if ~isempty(lr_node)
    fprintf('Bottom-right corner: Ux=%.4f mm, Uy=%.4f mm\n', u_x(lr_node)*1e3, u_y(lr_node)*1e3);
end

% Top-left corner
ul_node = find(abs(nodes(:,1))<1e-6 & abs(nodes(:,2)-H)<1e-6);
if ~isempty(ul_node)
    fprintf('Top-left corner: Ux=%.4f mm, Uy=%.4f mm\n', u_x(ul_node)*1e3, u_y(ul_node)*1e3);
end

% Center point
center_node = find(abs(nodes(:,1)-L/2)<1e-3 & abs(nodes(:,2)-H/2)<1e-3);
if ~isempty(center_node)
    fprintf('Center point: Ux=%.4f mm, Uy=%.4f mm\n', u_x(center_node)*1e3, u_y(center_node)*1e3);
end