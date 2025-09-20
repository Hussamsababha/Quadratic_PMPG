%% 
clc; clear; close all

%% fixed pars:
load constraintsFD.mat
nx = 200;
ny = 200; 
n = 201; 
tf = 5; 
dt = 0.0005; 
lx = 1; 
ly = 1; 
nsteps = 20; 
Re = 1e2; 
n_vars = n^2*2;
courant = 2.0;
%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;
x = linspace(0,lx,nx+1); hx = lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
[X,Y] = meshgrid(y,x);
h = lx/nx; 
%-----------------------------------------------------------------------
nu = 1/Re; 
U = zeros(n); V = zeros(n); 
C0 = 15; 
U(n, :) = 1 ;% - ((cosh(C0.*(X(end, :) - 0.5)))/(cosh(0.5*C0)));
x0 = [];

% plot(X(end, :), U(end, :))

% % Define video writer objects
% velocityVideo = VideoWriter('VelocityVideo_Re100.avi');  % Velocity video file
% accelerationVideo = VideoWriter('AccelerationVideo_Req100.avi');  % Acceleration video file
% 
% % Set frame rate (optional)
% velocityVideo.FrameRate = 10;
% accelerationVideo.FrameRate = 10;
% 
% % Open video files for writing
% open(velocityVideo);
% open(accelerationVideo);

% Define folder paths for saving figures
velocityFolder = 'velocityMidD_figures';  % Folder for velocity plots
accelerationFolder = 'accelerationMidD_figures';  % Folder for acceleration plots

% Create folders if they do not exist
if ~exist(velocityFolder, 'dir')
    mkdir(velocityFolder);
end
if ~exist(accelerationFolder, 'dir')
    mkdir(accelerationFolder);
end

% x0 = zeros(n_vars, 1);
for i = 1:nt
disp(['i = ' num2str(i)])
[Ux, Uy] = gradient(U, h);
[Vx, Vy] = gradient(V, h);
[Uxx, ~] = gradient(Ux, h);
[Vxx, ~] = gradient(Vx, h);
[~, Uyy] = gradient(Uy, h);
[~, Vyy] = gradient(Vy, h);
a = U.*Ux + V.*Uy - nu.*(Uxx + Uyy);
b = U.*Vx + V.*Vy - nu.*(Vxx + Vyy);
a = a';
b = b'; 
f = 2.*[a(:); b(:)];


%% Solving: 
options = optimoptions('quadprog','Display', 'none', 'Algorithm', 'interior-point-convex');
 %'TolFun', 1e-2, 'TolCon', 1e-2
tic
[optimal_Ut, fval, exitflag, output, lambda]  = quadprog(2*speye(n_vars), f, [], [], Aeq, beq, [], [], x0, options);
ans = toc
ut = optimal_Ut(1: n^2); vt = optimal_Ut(n^2 + 1: end);
Ut = reshape(ut, n, n)';
Vt = reshape(vt, n, n)';
x0 = optimal_Ut;
Ut_mag = sqrt(Ut.^2 + Vt.^2);
maxUt(i) = max(Ut_mag(:));
disp(['Max Ut = ' num2str(maxUt(i))])
dt(i) = (courant * hx/(maxUt(i)));
disp(['The time step = ' num2str(dt(i))])

%% compute the pressure gradient integral
momentum = (Ut + U.*Ux + V.*Uy - nu.*(Uxx + Uyy)).^2 + (Vt + U.*Vx + V.*Vy - nu.*(Vxx + Vyy)).^2; 
PMPG(i) = h^2*sum(momentum(:));
%% March in Time and compute a and b.
U = U + dt(i).* Ut;
V = V + dt(i).* Vt;
U_mag = sqrt(U.^2 + V.^2);


% if i==1|floor(nsteps*i/nt)>floor(nsteps*(i-1)/nt)
% %         % Create the velocity plot
%         figure;
%         surf(X, Y, U, 'LineStyle', 'none');
%         view(2);
%         % title(['Iteration ', num2str(i), ' - Velocity']);
%         colormap jet
% %         caxis([0 1])
%         % Save the figure as a .fig file
% %         savefig(fullfile(velocityFolder, ['velocity_iteration_', num2str(i), '.fig']));
% 
%         figure;
%         surf(X, Y, V, 'LineStyle', 'none');
%         view(2);
%         % title(['Iteration ', num2str(i), ' - Velocity']);
%         colormap jet
% %         caxis([0 1])
%         % Save the figure as a .fig file
% %         savefig(fullfile(velocityFolder, ['velocity_iteration_', num2str(i), '.fig']));
% 
        figure;
        surf(X, Y, sqrt(U.^2 + V.^2), 'LineStyle', 'none');
        view(2);
        title(['Iteration ', num2str(i), ' - Local Acceleration']);
        colormap jet
        colorbar
      % caxis([0 5])
% %         savefig(fullfile(accelerationFolder, ['acceleration_iteration_', num2str(i), '.fig']));
% %         % Capture the frame
% %         frame = getframe(gcf);
% %         writeVideo(velocityVideo, frame);  % Write to velocity video
% %         close(gcf);  % Close the figure to save memory
% %         % Create the line plot
% %         figure;
% %         plot(Y(:, 50), Ut_mag(:, 100));
% %         view(2);
% % %       title(['Iteration ', num2str(i), ' - Local Acceleration']);
% %         colormap jet
% % %       colorbar
% % %         caxis([0 5])
%         % Save the figure as a .fig file
% %         close all
% %         % Capture the frame
% %         frame = getframe(gcf);
% %         writeVideo(accelerationVideo, frame);  % Write to acceleration video
% %         close(gcf);  % Close the figure to save memory
% end

% if i==1|floor(nsteps*i/nt)>floor(nsteps*(i-1)/nt)
%     % Create surf plot
% %     clf
%     surf(X, Y, U_mag, 'EdgeColor', 'none'); % Surface plot of velocity magnitude
%     colormap jet
%     colorbar; % Add colorbar to indicate magnitude
%     axis equal; % Ensure equal scaling for axes
%     title(sprintf('Velocity Magnitude at t = %0.2g', i * dt)); % Add a title
%     xlabel('x'); ylabel('y'); zlabel('Velocity Magnitude'); % Label axes
%     view(2); % View from the top (2D)
%     drawnow; % Update plot
% end

end

Time = 0:0.001:1000; 
plot(Time(1:length(PMPG)), PMPG)
% ylim([5 15])
xlabel('Time (s)')
ylabel('$\mathcal{A}$', 'Interpreter','latex')
close(velocityVideo);
close(accelerationVideo);

%%
% Define starting points for particles
[startX, startY] = meshgrid(linspace(0,1,20), linspace(0,1,20)); % 20x20 grid

% Plot flow pathlines (steady flow)
figure;
streamline(X, Y, U, V, startX, startY)
axis equal
xlabel('X'); ylabel('Y')
title('Flow Pathlines')
grid on