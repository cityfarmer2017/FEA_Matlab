% Finite Element Analysis: Thin Plate (Plane Stress)
clear; clc;

%% 1. Mesh: Nodes & Elements
node_coords = [0, 0; 0.5, 0; 1.0, 0; 1.5, 0;...
               0, 0.5; 0.5, 0.5; 1.0, 0.5; 1.5, 0.5]; % 8 nodes
elements = [1,2,6,5; 2,3,7,6; 3,4,8,7]; % 3 elements (4-node quads)
num_elem = size(elements,1); num_nodes = size(node_coords,1);
dof_per_node = 2; total_dof = num_nodes*dof_per_node;

%% 2. Properties
E = 200e9; nu = 0.3; t = 0.02; p = 1500; % N/m (right-edge load)
D = (E/(1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2]; % Plane Stress D

%% 3. Boundary Conditions (Fixed DOFs: u=0/v=0)
fixed_dof = [1,2,4,6,8,9,10]; % Left edge (1,5: DOFs 1,2,9,10); Bottom edge (1-4: DOFs 2,4,6,8)

%% 4. Assemble Global Stiffness [K]
K = sparse(total_dof, total_dof);
gauss_pts = [-1/sqrt(3), 1/sqrt(3)]; gauss_wts = [1,1]; % 2x2 Gauss

for e = 1:num_elem
    elem_nodes = elements(e,:);
    x = node_coords(elem_nodes,1); y = node_coords(elem_nodes,2);
    elem_dof = [2*elem_nodes-1, 2*elem_nodes]; elem_dof = elem_dof(:);
    Ke = zeros(8,8);
    
    for i=1:2
        xi = gauss_pts(i); wxi = gauss_wts(i);
        for j=1:2
            eta = gauss_pts(j); wet = gauss_wts(j);
            
            % Shape functions/derivatives (4-node quad)
            N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
            dN_dxi = 0.25*[-(1-eta), (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25*[-(1-xi), -(1+xi), (1+xi), (1-xi)];
            
            % Jacobian
            J = [dN_dxi*x, dN_dxi*y; dN_deta*x, dN_deta*y];
            detJ = det(J); invJ = inv(J);
            
            % dN/dx, dN/dy
            dN_dx = invJ(1,1)*dN_dxi + invJ(1,2)*dN_deta;
            dN_dy = invJ(2,1)*dN_dxi + invJ(2,2)*dN_deta;
            
            % Strain matrix [B]
            B = zeros(3,8);
            for n=1:4
                B(1,2*n-1) = dN_dx(n);
                B(2,2*n) = dN_dy(n);
                B(3,2*n-1) = dN_dy(n); B(3,2*n) = dN_dx(n);
            end
            
            Ke = Ke + wxi*wet*detJ * (B'*D*B);
        end
    end
    Ke = Ke * t; % Multiply by thickness
    K(elem_dof,elem_dof) = K(elem_dof,elem_dof) + Ke;
end

%% 5. Assemble Load Vector {F}
F = zeros(total_dof,1);
% Right edge (nodes 4,8: DOFs 7 (u4), 15 (u8))
F(7) = 375; F(15) = 375; % Integral of p*0.5*(1-s)/s ds = 375 N each (total 750 N)

%% 6. Apply BCs & Solve
free_dof = setdiff(1:total_dof, fixed_dof);
u_free = K(free_dof,free_dof) \ F(free_dof);
u = zeros(total_dof,1); u(free_dof) = u_free;

%% 7. Print Results
fprintf('Node Displacements (u [m], v [m]):\n');
for n=1:num_nodes
    fprintf('Node %d: u = %.6e, v = %.6e\n', n, u(2*n-1), u(2*n));
end