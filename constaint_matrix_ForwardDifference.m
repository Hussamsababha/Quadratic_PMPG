clc; clear

% Grid size
n = 200; % nxn grid
n_points = n^2; % Total grid points
n_vars = 2 * n_points; % u and v variables for each point
% Problem parameters
dx = 1; % Grid spacing in x
dy = 1; % Grid spacing in y

%%
% divergence constraint
A_cont = sparse(0, n_vars); % sparse representation 
b_cont = []; % Right-hand side (zeros)

% Loop over interior points (excluding boundaries)
for i = 1:n-1
    for j = 1:n-1
        row = zeros(1, n_vars); % Initialize row for this equation
        % mid point difference
        row((i-1)*n + (j + 1)) = 1 / dx; % u(i+1,j)
        row((i-1)*n + (j )) = -1 / dx; % u(i,j)
        row(n_points + (i-1)*n + j) = -1 / dy; % v(i,j+1)
        row(n_points + (i)*n + j) = +1 / dy;   % v(i,j)
        
        % Add this row to A_cont
        A_cont = [A_cont; sparse(row)]; % sparse representation 
        b_cont = [b_cont; 0]; % RHS is zero
    end
end

disp(['[Size of A_cont]' num2str(size(A_cont))])
% Boundary constraints
A_bc_u = sparse(0, n_vars); % Initialize as zero matrix
A_bc_v = sparse(0, n_vars); % Initialize as zero matrix

% Loop over boundary points
for i = 1:n
    for j = 1:n
        if i == 1 || i == n || j == 1 || j == n % Boundary points
            % u = 0 (boundary condition for u)
            row_u = zeros(1, n_vars);
            row_u((i-1)*n + j) = 1; % u(i,j)
            A_bc_u = [A_bc_u; sparse(row_u)];
            % v = 0 (boundary condition for v)
            row_v = zeros(1, n_vars);
            row_v(n_points + (i-1)*n + j) = 1; % v(i,j)
            %A_bc = [A_bc; row_v];
            A_bc_v = [A_bc_v; sparse(row_v)];
        end
    end
end

disp(['[Size of A_bc_u]' num2str(size(A_bc_u))])
disp(['[Size of A_bc_v]' num2str(size(A_bc_v))])

A_bc = [A_bc_u; A_bc_v];
Aeq = [A_cont; A_bc]; % Combine the constraints
beq = zeros(size(Aeq, 1), 1); 
beq = sparse(beq);

save('constraintsFD', 'Aeq', 'beq')
% A = full(Aeq); 
% b = full(beq); 
% 
% %% Objective function
% % M = 2 * eye(n_vars); % Diagonal matrix for u and v terms
% % c = 2 * [a(:); b(:)]; % Vector of coefficients for linear terms
% 
% syms_list_u = []; % Initialize for u variables
% syms_list_v = []; % Initialize for v variables
% 
% 
% % Generate u variables
% for i = 1:n
%     for j = 1:n
%         syms(sprintf('u%d%d', i, j)); % Dynamically create symbolic variables
%         syms_list_u = [syms_list_u; sym(sprintf('u%d%d', i, j))];
%     end
% end
% 
% % Generate v variables
% for i = 1:n
%     for j = 1:n
%         syms(sprintf('v%d%d', i, j)); % Dynamically create symbolic variables
%         syms_list_v = [syms_list_v; sym(sprintf('v%d%d', i, j))];
%     end
% end
% 
% % Combine into one vector if needed
% syms_vector = [syms_list_u; syms_list_v];
% 
% constraints = A * syms_vector
% 
