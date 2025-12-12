%% 薄板变形平面应力问题有限元分析
% 开发人员：[您的姓名/小组名称]
% 日期：[当前日期]

clear all; close all; clc;

%% 1. 问题参数定义
E = 200e9;          % 弹性模量 (Pa)
nu = 0.3;           % 泊松比
t = 0.02;           % 板厚度 (m)
L = 1.5;            % 板长度 (m)
H = 0.5;            % 板高度 (m)
p = 1500;           % 分布载荷 (N/m)

% 转换为面力 (Pa)
traction_y = -p/t;  % 负号表示向下

%% 2. 网格生成
nx = 30;            % x方向单元数
ny = 10;            % y方向单元数
nnodes = (nx+1)*(ny+1);     % 节点总数
nelems = nx*ny;             % 单元总数

% 生成节点坐标
x = linspace(0, L, nx+1);
y = linspace(0, H, ny+1);
[X, Y] = meshgrid(x, y);
nodes = [X(:), Y(:)];

% 生成单元连接性（4节点四边形单元）
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

%% 3. 材料矩阵（平面应力）
D = (E/(1-nu^2)) * [1, nu, 0;
                   nu, 1, 0;
                   0, 0, (1-nu)/2];

%% 4. 全局刚度矩阵和力向量组装
ndof = 2 * nnodes;          % 总自由度
K = sparse(ndof, ndof);     % 稀疏存储
F = zeros(ndof, 1);

% 高斯积分点（2×2）
gauss_points = [-1/sqrt(3), -1/sqrt(3);
                 1/sqrt(3), -1/sqrt(3);
                 1/sqrt(3),  1/sqrt(3);
                -1/sqrt(3),  1/sqrt(3)];
weights = [1, 1, 1, 1];

% 单元循环
for e = 1:nelems
    % 获取节点坐标
    node_ids = elements(e, :);
    elem_nodes = nodes(node_ids, :);
    
    % 初始化单元刚度矩阵
    ke = zeros(8, 8);
    
    % 数值积分
    for gp = 1:4
        xi = gauss_points(gp, 1);
        eta = gauss_points(gp, 2);
        
        % 形函数及其导数
        N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), ...
                    (1+xi)*(1+eta), (1-xi)*(1+eta)];
        
        dN_dxi = 0.25 * [-(1-eta),  (1-eta),  (1+eta), -(1+eta)];
        dN_deta = 0.25 * [-(1-xi), -(1+xi),  (1+xi),  (1-xi)];
        
        % Jacobian矩阵
        J = [dN_dxi * elem_nodes(:,1), dN_dxi * elem_nodes(:,2);
             dN_deta * elem_nodes(:,1), dN_deta * elem_nodes(:,2)];
        
        detJ = det(J);
        invJ = inv(J);
        
        % 应变-位移矩阵B
        B = zeros(3, 8);
        for i = 1:4
            dN_dx = invJ(1,1)*dN_dxi(i) + invJ(1,2)*dN_deta(i);
            dN_dy = invJ(2,1)*dN_dxi(i) + invJ(2,2)*dN_deta(i);
            
            B(1, 2*i-1) = dN_dx;
            B(2, 2*i) = dN_dy;
            B(3, 2*i-1) = dN_dy;
            B(3, 2*i) = dN_dx;
        end
        
        % 单元刚度矩阵贡献
        ke = ke + B' * D * B * detJ * weights(gp) * t;
    end
    
    % 组装到全局矩阵
    dof_indices = zeros(8, 1);
    for i = 1:4
        dof_indices(2*i-1:2*i) = [2*node_ids(i)-1, 2*node_ids(i)];
    end
    
    K(dof_indices, dof_indices) = K(dof_indices, dof_indices) + ke;
end

%% 5. 施加载荷（顶部边缘）
% 找出顶部边缘的节点
top_nodes = find(nodes(:,2) == H);

% 计算节点间距以分配载荷
for i = 1:length(top_nodes)
    node_id = top_nodes(i);
    
    % 确定相邻节点
    if i == 1
        dx = nodes(top_nodes(2),1) - nodes(node_id,1);
    elseif i == length(top_nodes)
        dx = nodes(node_id,1) - nodes(top_nodes(end-1),1);
    else
        dx = (nodes(top_nodes(i+1),1) - nodes(top_nodes(i-1),1))/2;
    end
    
    % 等效节点力（线性分布）
    F(2*node_id) = F(2*node_id) + traction_y * dx * t / 2;
end

%% 6. 边界条件（简支边：左右两侧）
% 找出左右边缘的节点
left_nodes = find(nodes(:,1) == 0);
right_nodes = find(nodes(:,1) == L);

% 约束自由度：左右边缘u=0, v自由
constrained_dofs = [];
for i = 1:length(left_nodes)
    constrained_dofs = [constrained_dofs, 2*left_nodes(i)-1];
end
for i = 1:length(right_nodes)
    constrained_dofs = [constrained_dofs, 2*right_nodes(i)-1];
end

constrained_dofs = unique(constrained_dofs);

%% 7. 应用边界条件并求解
free_dofs = setdiff(1:ndof, constrained_dofs);

% 缩减系统求解
U = zeros(ndof, 1);
U(free_dofs) = K(free_dofs, free_dofs) \ F(free_dofs);

% 提取位移
u_x = U(1:2:end);   % x方向位移
u_y = U(2:2:end);   % y方向位移

%% 8. 结果分析和可视化
fprintf('=== 有限元分析结果 ===\n');
fprintf('节点总数: %d\n', nnodes);
fprintf('单元总数: %d\n', nelems);
fprintf('最大垂直位移: %.6e m\n', min(u_y));
fprintf('最大水平位移: %.6e m\n', max(abs(u_x)));

% 创建结果表格
results_table = table((1:nnodes)', nodes(:,1), nodes(:,2), ...
                     u_x, u_y, ...
                     'VariableNames', {'Node', 'X', 'Y', 'Ux', 'Uy'});

% 输出前20个节点的结果
fprintf('\n前20个节点的位移结果:\n');
disp(results_table(1:min(20, nnodes), :));

% 沿中心线(y=H/2)的位移
centerline_nodes = find(abs(nodes(:,2) - H/2) < 1e-6);
[~, idx] = sort(nodes(centerline_nodes, 1));
centerline_nodes = centerline_nodes(idx);

% 绘图
figure('Position', [100, 100, 1200, 500]);

% 子图1: x方向位移
subplot(1, 2, 1);
plot(nodes(centerline_nodes, 1), u_x(centerline_nodes)*1e3, 'b-o', 'LineWidth', 2);
xlabel('x坐标 (m)');
ylabel('水平位移 U_x (mm)');
title('沿中心线的水平位移分布');
grid on;
xlim([0, L]);

% 子图2: y方向位移
subplot(1, 2, 2);
plot(nodes(centerline_nodes, 1), u_y(centerline_nodes)*1e3, 'r-s', 'LineWidth', 2);
xlabel('x坐标 (m)');
ylabel('垂直位移 U_y (mm)');
title('沿中心线的垂直位移分布');
grid on;
xlim([0, L]);

% 保存结果到文件
writetable(results_table, 'displacement_results.csv');
save('FEM_results.mat', 'nodes', 'u_x', 'u_y', 'elements');

fprintf('\n结果已保存到文件:\n');
fprintf('  - displacement_results.csv (位移数据表)\n');
fprintf('  - FEM_results.mat (MATLAB数据文件)\n');

%% 9. 附加分析：位移等高线图
figure('Position', [100, 100, 1000, 400]);

% 创建网格用于等高线
x_grid = reshape(nodes(:,1), [ny+1, nx+1]);
y_grid = reshape(nodes(:,2), [ny+1, nx+1]);
u_x_grid = reshape(u_x, [ny+1, nx+1]);
u_y_grid = reshape(u_y, [ny+1, nx+1]);

subplot(1, 2, 1);
contourf(x_grid, y_grid, u_x_grid*1e3, 20);
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('水平位移 U_x (mm)');
axis equal;

subplot(1, 2, 2);
contourf(x_grid, y_grid, u_y_grid*1e3, 20);
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('垂直位移 U_y (mm)');
axis equal;