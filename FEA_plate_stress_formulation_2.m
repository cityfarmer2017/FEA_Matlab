%% 薄板平面应力分析 - 左侧和底边简支，水平分布载荷
% 开发人员：[您的姓名/小组名称]
% 日期：[当前日期]

clear all; close all; clc;

%% 1. 问题参数定义
E = 200e9;          % 弹性模量 (Pa)
nu = 0.3;           % 泊松比
t = 0.02;           % 板厚度 (m)
L = 1.5;            % 板长度 (m)
H = 0.5;            % 板高度 (m)
p = 1500;           % 分布载荷 (N/m) - 水平方向从左到右

% 转换为线载荷强度 (N/m)
traction_x = p;     % 水平向右的正方向

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

%% 5. 施加载荷（右侧边缘 - 水平方向）
% 找出右侧边缘的节点（x = L）
right_nodes = find(abs(nodes(:,1) - L) < 1e-6);

% 按y坐标排序
[~, idx] = sort(nodes(right_nodes, 2));
right_nodes = right_nodes(idx);

% 计算节点间距以分配载荷
for i = 1:length(right_nodes)
    node_id = right_nodes(i);
    
    % 确定相邻节点（用于计算载荷分配长度）
    if i == 1
        % 第一个节点，只考虑上侧一半
        dy = nodes(right_nodes(2),2) - nodes(node_id,2);
        dy = dy / 2;
    elseif i == length(right_nodes)
        % 最后一个节点，只考虑下侧一半
        dy = nodes(node_id,2) - nodes(right_nodes(end-1),2);
        dy = dy / 2;
    else
        % 中间节点，考虑两侧各一半
        dy = (nodes(right_nodes(i+1),2) - nodes(right_nodes(i-1),2)) / 2;
    end
    
    % 等效节点力（水平向右）
    F(2*node_id-1) = F(2*node_id-1) + traction_x * dy;
end

%% 6. 边界条件（简支边：左侧和底边）
% 找出左侧边缘的节点（x = 0）
left_nodes = find(abs(nodes(:,1)) < 1e-6);

% 找出底边边缘的节点（y = 0）
bottom_nodes = find(abs(nodes(:,2)) < 1e-6);

% 约束自由度：
% 左侧边：约束水平位移 (u=0)
% 底边：约束垂直位移 (v=0)
constrained_dofs = [];

% 左侧边节点：约束x方向位移（自由度序号为奇数）
for i = 1:length(left_nodes)
    constrained_dofs = [constrained_dofs, 2*left_nodes(i)-1];
end

% 底边节点：约束y方向位移（自由度序号为偶数）
for i = 1:length(bottom_nodes)
    constrained_dofs = [constrained_dofs, 2*bottom_nodes(i)];
end

% 去除重复（左下角节点同时属于左侧和底边）
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
fprintf('最大水平位移: %.6e m\n', max(u_x));
fprintf('最小水平位移: %.6e m\n', min(u_x));
fprintf('最大垂直位移: %.6e m\n', max(u_y));
fprintf('最小垂直位移: %.6e m\n', min(u_y));
fprintf('右下角水平位移: %.6e m\n', u_x(right_nodes(end)));
fprintf('右上角水平位移: %.6e m\n', u_x(right_nodes(1)));

% 创建结果表格
results_table = table((1:nnodes)', nodes(:,1), nodes(:,2), ...
                     u_x, u_y, ...
                     'VariableNames', {'Node', 'X', 'Y', 'Ux', 'Uy'});

% 输出前20个节点的结果
fprintf('\n前20个节点的位移结果:\n');
disp(results_table(1:min(20, nnodes), :));

% 计算沿特定路径的位移
% 沿顶边 (y=H) 的位移
top_edge_nodes = find(abs(nodes(:,2) - H) < 1e-6);
[~, idx] = sort(nodes(top_edge_nodes, 1));
top_edge_nodes = top_edge_nodes(idx);

% 沿右侧边 (x=L) 的位移
[~, idx] = sort(nodes(right_nodes, 2));
right_edge_nodes = right_nodes(idx);

%% 9. 可视化结果
figure('Position', [50, 50, 1400, 900]);

% 子图1: 沿顶边的位移
subplot(2, 3, 1);
plot(nodes(top_edge_nodes, 1), u_x(top_edge_nodes)*1e3, 'b-o', 'LineWidth', 2, 'DisplayName', 'U_x');
hold on;
plot(nodes(top_edge_nodes, 1), u_y(top_edge_nodes)*1e3, 'r-s', 'LineWidth', 2, 'DisplayName', 'U_y');
xlabel('x坐标 (m)');
ylabel('位移 (mm)');
title('沿顶边 (y=H) 的位移分布');
legend('位置', 'best');
grid on;
xlim([0, L]);

% 子图2: 沿右侧边的位移
subplot(2, 3, 2);
plot(nodes(right_edge_nodes, 2), u_x(right_edge_nodes)*1e3, 'b-o', 'LineWidth', 2, 'DisplayName', 'U_x');
hold on;
plot(nodes(right_edge_nodes, 2), u_y(right_edge_nodes)*1e3, 'r-s', 'LineWidth', 2, 'DisplayName', 'U_y');
xlabel('y坐标 (m)');
ylabel('位移 (mm)');
title('沿右侧边 (x=L) 的位移分布');
legend('位置', 'best');
grid on;
xlim([0, H]);

% 子图3: 水平位移云图
subplot(2, 3, 3);
% 重新整形为网格格式
x_grid = reshape(nodes(:,1), [ny+1, nx+1]);
y_grid = reshape(nodes(:,2), [ny+1, nx+1]);
u_x_grid = reshape(u_x, [ny+1, nx+1]);
contourf(x_grid, y_grid, u_x_grid*1e3, 20, 'LineStyle', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('水平位移 U_x (mm)');
axis equal tight;

% 子图4: 垂直位移云图
subplot(2, 3, 4);
u_y_grid = reshape(u_y, [ny+1, nx+1]);
contourf(x_grid, y_grid, u_y_grid*1e3, 20, 'LineStyle', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('垂直位移 U_y (mm)');
axis equal tight;

% 子图5: 位移矢量图（变形放大）
subplot(2, 3, 5);
scale_factor = 0.1 / max(sqrt(u_x.^2 + u_y.^2)); % 适当缩放
quiver(nodes(:,1), nodes(:,2), u_x*scale_factor, u_y*scale_factor, 0, 'b');
hold on;
plot(nodes(:,1), nodes(:,2), 'k.', 'MarkerSize', 8);
xlabel('x (m)');
ylabel('y (m)');
title('位移矢量图（缩放显示）');
axis equal tight;
grid on;

% 子图6: 变形前后的网格对比
subplot(2, 3, 6);
% 原始网格
for e = 1:nelems
    node_ids = elements(e, :);
    patch(nodes(node_ids, 1), nodes(node_ids, 2), 'w', 'EdgeColor', [0.7 0.7 0.7], 'FaceColor', 'none');
end
hold on;
% 变形后网格（放大显示）
deformation_scale = 100; % 放大倍数以便观察
deformed_nodes = nodes + [u_x, u_y] * deformation_scale;
for e = 1:min(nelems, 200) % 只绘制部分单元以免太密
    node_ids = elements(e, :);
    patch(deformed_nodes(node_ids, 1), deformed_nodes(node_ids, 2), 'w', 'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 1.5);
end
xlabel('x (m)');
ylabel('y (m)');
title(sprintf('变形前后网格对比（变形放大%d倍）', deformation_scale));
legend('原始网格', '变形网格', '位置', 'best');
axis equal tight;
grid on;

% 保存结果到文件
writetable(results_table, 'displacement_results_horizontal_load.csv');
save('FEM_results_horizontal_load.mat', 'nodes', 'u_x', 'u_y', 'elements', 'right_nodes', 'top_edge_nodes');

fprintf('\n结果已保存到文件:\n');
fprintf('  - displacement_results_horizontal_load.csv (位移数据表)\n');
fprintf('  - FEM_results_horizontal_load.mat (MATLAB数据文件)\n');

%% 10. 附加统计信息
fprintf('\n=== 位移统计信息 ===\n');
fprintf('水平位移范围: %.4f 到 %.4f mm\n', min(u_x)*1e3, max(u_x)*1e3);
fprintf('垂直位移范围: %.4f 到 %.4f mm\n', min(u_y)*1e3, max(u_y)*1e3);
fprintf('平均水平位移: %.4f mm\n', mean(u_x)*1e3);
fprintf('平均垂直位移: %.4f mm\n', mean(u_y)*1e3);

% 计算特定点的位移
% 右上角
ur_node = find(abs(nodes(:,1)-L)<1e-6 & abs(nodes(:,2)-H)<1e-6);
if ~isempty(ur_node)
    fprintf('右上角位移: Ux=%.4f mm, Uy=%.4f mm\n', u_x(ur_node)*1e3, u_y(ur_node)*1e3);
end

% 右下角
lr_node = find(abs(nodes(:,1)-L)<1e-6 & abs(nodes(:,2))<1e-6);
if ~isempty(lr_node)
    fprintf('右下角位移: Ux=%.4f mm, Uy=%.4f mm\n', u_x(lr_node)*1e3, u_y(lr_node)*1e3);
end

% 左上角
ul_node = find(abs(nodes(:,1))<1e-6 & abs(nodes(:,2)-H)<1e-6);
if ~isempty(ul_node)
    fprintf('左上角位移: Ux=%.4f mm, Uy=%.4f mm\n', u_x(ul_node)*1e3, u_y(ul_node)*1e3);
end