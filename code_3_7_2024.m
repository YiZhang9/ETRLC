clc
close all
clear

T = 2498; eta = 0.1; mu = 1; lambda = 1.2; epsilon = 0.41; rho = 0.35; sigma = 1e-5;
gamma0 = 1; varphi = 0.5; v = 20;
m = 3; n = 2; varpi = -0.1; K = 0.9; omega = -2;
boundary = 0.001; line_value = 3; N = 8;
a = [2;2;3;5;2;2;3;5]; b = [1;2;3;2;2;4;2;3];
W = diag([1 1 -1 1 1 1 -1 1]);
S = diag([m m n m m m n m]);


n_a = 40; max_duration = 50;
DoS = ones(1, T); DoS_durations = zeros(1, n_a);
for i = 1:n_a
    DoS_durations(i) = randi([1, max_duration]);
    DoSstartIdx = randi([1, T - DoS_durations(i) + 1], 1, 1);
    while any(DoS(DoSstartIdx:DoSstartIdx + DoS_durations(i) - 1) == 0)
        DoSstartIdx = randi([1, T - DoS_durations(i) + 1], 1, 1);
    end
    DoS(DoSstartIdx:DoSstartIdx + DoS_durations(i) - 1) = 0;
end
Xi_a = sum(DoS_durations);


FDI = zeros(N, T);
for k = 1:T
    for i = 1:4
        FDI(i, k) = 2 * rand() - 1;
    end
    for i = 5:8
        FDI(i, k) = sin(i*k);
    end
end

y1 = zeros(N, T); u1 = zeros(N, T); phi_hat1 = zeros(N, T);
Delta_u1 = zeros(N, T); Delta_y1 = zeros(N, T);
xi1 = zeros(N, T); e_y1 = zeros(N, T); Delta_xi1 = zeros(N, T);
y1(:, 1) = 0.5; phi_hat1(:, :, 1) = 1;

y_d = zeros(1, T);




for k =1:1400
    y_d(k) = 6;
end
for k =1400:T
    y_d(k) = 7*sin(pi*k/2200) - 8*cos(pi*k/2400);
end



k = 1;
for i=1:N
    y1(i, k+1) = y1(i, k)*u1(i, k)/(1+y1(i, k)^a(i))+b(i)*u1(i, k);
end
e_y1(:, k) = ones(N, 1)*y_d(k)-W*S\y1(:, k);
xi1(1, k) = y_d(1, k) - 1/m*y1(1, k) - y1(2, k) - 1/m*n*y1(1, k);
xi1(2, k) = - y1(4, k) - 1/n*m*y1(2, k);
xi1(3, k) = - y1(1, k) - 1/n*m*y1(3, k);
xi1(4, k) = - y1(6, k) - 1/m*n*y1(4, k);
xi1(5, k) = y_d(1, k) - 1/m*y1(5, k);
xi1(6, k) =  - y1(5, k) - 1/n*m*y1(6, k) - y_d(1, k) - 1/m*y1(6, k);
xi1(7, k) = - y1(6, k) - 1/m*n*y1(7, k);
xi1(8, k) = y1(7, k) - y1(8, k);
Delta_xi11(:, k) = xi1(:, k);

for k = 2:T-1
    xi1(1, k) = y_d(1, k) - 1/m*y1(1, k) - y1(2, k) - 1/m*n*y1(1, k);
    xi1(2, k) = - y1(4, k) - 1/n*m*y1(2, k);
    xi1(3, k) = - y1(1, k) - 1/n*m*y1(3, k);
    xi1(4, k) = - y1(6, k) - 1/m*n*y1(4, k);
    xi1(5, k) = y_d(1, k) - 1/m*y1(5, k);
    xi1(6, k) =  - y1(5, k) - 1/n*m*y1(6, k) - y_d(1, k) - 1/m*y1(6, k);
    xi1(7, k) = - y1(6, k) - 1/m*n*y1(7, k);
    xi1(8, k) = y1(7, k) - y1(8, k);
    for i = 1:N
        phi_hat1(i, k) = phi_hat1(i, k-1) ...
            + eta * (xi1(i, k)-xi1(i, k-1)) * Delta_u1(i, k-1) ...
            /(Delta_u1(i, k-1)^2+mu) ...
            - eta * phi_hat1(i, k-1)*Delta_u1(i, k-1)*Delta_u1(i, k-1)' ...
            /(Delta_u1(i, k-1)^2+mu);
        if abs(phi_hat1(i, k)) < sigma ...
                || abs(Delta_u1(i, k)) < sigma ...
                || sign(phi_hat1(i, k)) ~= sign(phi_hat1(i, 1))
            phi_hat1(i, k) = phi_hat1(i, 1);
        end
        Delta_u1(i, k) = eta*phi_hat1(i, k)' ...
            / (varpi+phi_hat1(i, k)^2) ...
            *xi1(i, k);
        u1(i, k) = u1(i, k-1) + Delta_u1(i, k);
    end
    for i=1:N
        y1(i, k+1) = y1(i, k)*u1(i, k)/(1+y1(i, k)^a(i))+b(i)*u1(i, k);
    end
    Delta_y1(:, k+1) = y1(i, k+1) - y1(i, k);
    e_y1(:, k) = ones(N, 1)*y_d(k)-W*S\y1(:, k);
end
varPsi = round(xi1/e_y1);

A = [0,1,-1,0,0,0,0,0;
    0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,1,0,0,0;
    0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,1;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,-1,0];
G = diag([1,0,0,0,0,1,0,0]);
% L_absG = S\W*(varPsi - G);
% L_absA = S*L_absG\S;
% D_hat = diag((W*S\A*S*W)*ones(N,1));
% L_absA2 = D_hat-W*A*W;
% L_absG2 = S\L_absA2*S;
% varPsi2 = W*S*L_absG2 + G;
% L_bar = W*L_absA*W;
% L_bar2 =


% W = diag([1,1,-1,1,1,1,-1,1]);
% S = diag([3,3,4,3,3,3,4,3]);
% A = [0,1,-1,0,0,0,0,0;
%      0,0,0,1,0,0,0,0;
%      0,0,0,0,0,0,0,0;
%      0,0,0,0,1,0,0,0;
%      0,0,0,0,0,1,0,0;
%      0,0,0,0,0,0,0,1;
%      0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,-1,0];
% absA = W*A*W;
% D_absA = diag(absA*ones(8, 1));
% L_absA = D_absA - absA;
% D_hat = diag((S\A*S)*ones(8, 1));
% L_absG = D_hat - S\A*S;
% G = diag([1,0,0,0,0,1,0,0]);
% varPsi = W*S*L_absG+G;


trigger = zeros(N, T);
y = zeros(N, T); u = zeros(N, T); phi_hat = zeros(N, T);
Delta_u = zeros(N, T); Delta_y = zeros(N, T); r = zeros(N, T);
xi_hat = zeros(N, T); xi = zeros(N, T); Delta_xi = zeros(N, T);
xi_tilde = zeros(N, T); theta = zeros(1, T); ki = ones(N, T);
y(:, 1) = 0.5; phi_hat(:, 1) = 1;
P = zeros(N, N, T); Phi = zeros(N, N, T); H = zeros(N, N, T);


k = 1;
for i=1:N
    y(i, k+1) = y(i, k)*u(i, k)/(1+y(i, k)^a(i))+b(i)*u(i, k);
end
xi(1, k) = y_d(1, k) - 1/m*y1(1, k) - y1(2, k) - 1/m*n*y1(1, k);
xi(2, k) = - y1(4, k) - 1/n*m*y1(2, k);
xi(3, k) = - y1(1, k) - 1/n*m*y1(3, k);
xi(4, k) = - y1(6, k) - 1/m*n*y1(4, k);
xi(5, k) = y_d(1, k) - 1/m*y1(5, k);
xi(6, k) =  - y1(5, k) - 1/n*m*y1(6, k) - y_d(1, k) - 1/m*y1(6, k);
xi(7, k) = - y1(6, k) - 1/m*n*y1(7, k);
xi(8, k) = y1(7, k) - y1(8, k);


% xi(:, k) = varPsi*(ones(N, 1)*y_d(k)-W*S\y(:, k));
y_ki = [y(1, ki(1)); y(2, ki(2)); y(3, ki(3)); y(4, ki(4)); ...
    y(5, ki(5)); y(6, ki(6)); y(7, ki(7)); y(8, ki(8))];
xi_tilde(:, k) = varPsi*(ones(N, 1)*y_d(k)-W*S\y_ki);
for i = 1:N
    H(i, i, k) = eta*phi_hat(i, k)' / (varpi+phi_hat(i, k)^2);
end
Phi(:, :, k) = diag(phi_hat(:, k));
Delta(:, k) = S\W*(y(:, k) - y_ki);
Delta_xi(:, k) = xi(:, k);


for k = 2:T-1
    xi(1, k) = y_d(1, k) - 1/m*y1(1, k) - y1(2, k) - 1/m*n*y1(1, k);
    xi(2, k) = - y1(4, k) - 1/n*m*y1(2, k);
    xi(3, k) = - y1(1, k) - 1/n*m*y1(3, k);
    xi(4, k) = - y1(6, k) - 1/m*n*y1(4, k);
    xi(5, k) = y_d(1, k) - 1/m*y1(5, k);
    xi(6, k) =  - y1(5, k) - 1/n*m*y1(6, k) - y_d(1, k) - 1/m*y1(6, k);
    xi(7, k) = - y1(6, k) - 1/m*n*y1(7, k);
    xi(8, k) = y1(7, k) - y1(8, k);
    xi(:, k) = DoS(1, k)*xi(:, k);
    y_ki = [y(1, ki(1)); y(2, ki(2)); y(3, ki(3)); y(4, ki(4)); ...
        y(5, ki(5)); y(6, ki(6)); y(7, ki(7)); y(8, ki(8))];
    xi_tilde(:, k) = varPsi*(ones(N, 1)*y_d(k)-W*S\y_ki);
    Delta(:, k) = S\W*(y(:, k) - y_ki);

    P(:, :, k-1) = - Phi(:, :, k-1) * H(:, :, k-1);
    min_singular_P = svds(P(:, :, k-1), 1, 'smallest');
    theta(1, k - 1) = norm(P(:, :, k-1))^2 * norm(varPsi)^2 ...
        /(min_singular_P*(min_singular_P^2 - min_singular_P * norm(P(:, :, k-1))^2));

    for i = 1:N
        r(i, k-1) = abs(xi_tilde(i, k-1)) - theta(1, k - 1) * abs(Delta(i, k-1));
        if  r(i, k-1) <= omega
            trigger(i, k) = 1;
            phi_hat(i, k) = phi_hat(i, k-1) ...
                + eta * (xi(i, k)-xi(i, k-1)) * Delta_u(i, k-1) ...
                /(Delta_u(i, k-1)^2+mu) ...
                - eta * phi_hat(i, ki(i))*Delta_u(i, k-1)*Delta_u(i, k-1)' ...
                /(Delta_u(i, k-1)^2+mu);
            if abs(phi_hat(i, k)) < sigma ...
                    || abs(Delta_u(i, k)) < sigma ...
                    || sign(phi_hat(i, k)) ~= sign(phi_hat(i, 1))
                phi_hat(i, k) = phi_hat(i, 1);
            end
            Delta_u(i, k) = DoS(1, k)*eta*phi_hat(i, k)' ...
                / (varpi+phi_hat(i, k)^2) ...
                *xi(i, k);
            u(i, k) = u(i, k-1) + Delta_u(i, k);
            ki(i, k) = k;
        else
            phi_hat(i, k) = phi_hat(i, k-1);
            u(i, k) = u(i, k-1);
            ki(i, k) = ki(i, k - 1);
        end
    end
    Phi(:, :, k) = diag(phi_hat(:, k));
    for i = 1:N
        H(i, i, k) = eta*phi_hat(i, k)'/ (varpi+phi_hat(i, k)^2);
    end
    for i=1:N
        y(i, k+1) = y(i, k)*u(i, k)/(1+y(i, k)^a(i))+b(i)*u(i, k) + FDI(i, k);
    end
    Delta_y(:, k+1) = y(i, k+1) - y(i, k);
end

markers = ["none";"none";"none";"none";"none";"none";".";"none";"*"];

colors2 = [
    0.0, 0.45, 0.74;  % 蓝色
    0.85, 0.33, 0.1;  % 红色
    0.93, 0.69, 0.13; % 黄色
    0.49, 0.18, 0.56; % 紫色
    0.47, 0.67, 0.19; % 绿色
    0.3, 0.75, 0.93;  % 天蓝色
    0.64, 0.08, 0.18; % 葡萄酒红
    0.5, 0.5, 0.5;];

% color2

colors1 = [0.0, 0.6, 0.9;    % 亮蓝
    1.0, 0.4, 0.0;    % 亮橙
    0.9, 0.9, 0.0;    % 亮黄
    0.6, 0.0, 0.6;    % 亮紫
    0.0, 0.8, 0.0;    % 亮绿
    0.0, 0.75, 1.0;   % 天蓝
    1.0, 0.0, 0.0;    % 亮红
    0.75, 0.75, 0.75];% 亮灰
linestyle = ["-","--","-.",":","-","--",":","-."];




d = 2;

fig1 = figure(1);
hold on
ylim([-d d]);
yRange = ylim;
for i = 1:T-1
    if DoS(i) == 0
        color = [0.95 0.95 0.95]; % 灰色
        patch([i, i+1, i+1, i], ...
            [yRange(1), yRange(1), yRange(2), yRange(2)], ...
            color, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end
hold on
for i = 1:4
    plot(FDI(i, :), 'Color', colors1(i,:),...
        'DisplayName', ['noise', num2str(i)]);
end
for i = 5:8
    plot(FDI(i, :), 'Color', colors1(i,:),...
        'DisplayName', ['FDI', num2str(i-4)]);
end
% legend('Interpreter', 'latex', 'Location', 'southeast');
xlabel('Time(k)', 'Interpreter', 'latex');
ylabel('Hybrid Cyber Attacks for different Agent', 'Interpreter', 'latex');
hold off;
set(gca, 'Layer', 'top');

sub_axes1 = axes('Position',[0.6 0.6 .25 .25]);
box on
hold on
for i = 1:T-1
    if DoS(i) == 0
        color = [0.95 0.95 0.95]; % 灰色
        patch([i, i+1, i+1, i], ...
            [yRange(1), yRange(1), yRange(2), yRange(2)], ...
            color, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

for i = 1:4
    plot(FDI(i, :), 'Color', colors1(i,:),...
        'DisplayName', ['noise', num2str(i)]);
end
for i = 5:8
    plot(FDI(i, :), 'Color', colors1(i,:),...
        'DisplayName', ['FDI', num2str(i-4)]);
end
xlim(sub_axes1, [500 600]);
ylim(sub_axes1, [-1 1]);


d = 70;

fig2 = figure(2);


% subplot_positions1 = [0.05, 0.55, 0.4, 0.4]; % 第一个小图的位置和大小
% subplot_positions2 = [0.55, 0.55, 0.4, 0.4]; % 第二个小图的位置和大小
% subplot_positions3 = [0.05, 0.1, 0.4, 0.4]; % 第三个小图的位置和大小
% subplot_positions4 = [0.55, 0.1, 0.4, 0.4];


positions1 = [0.25, 0.62, 0.01, 0.01];
% positions2 = [0.81, 0.755, 0.01, 0.01];
% positions3 = [0.39, 0.29, 0.01, 0.01];
positions4 = [0.9, 0.55, 0.01, 0.01];

subplot(2,2,1);
% sub1.Position = subplot_positions1;

plot(y_d, 'LineWidth', 1.5,...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '\boldmath$y_d$')
hold on
plot(m*y_d, 'LineWidth', 1.5, ...
    'DisplayName', '\boldmath$m y_d$');
plot(-n*y_d, 'LineWidth', 1.5,...
    'DisplayName', '\boldmath$-ny_d$');
xlim([1, 2500])
ylim([-d d]);
% xlabel('Time(k)', 'Interpreter', 'latex');
% ylabel('Asymmetric Consensus Value', 'Interpreter', 'latex');
lgd1=legend('Interpreter', 'latex', 'Location', positions1);
lgd1.FontSize =6;
title('(a)');
hold off




subplot(2,2,3);
title('(c)');
% sub2.Position = subplot_positions3;

hold on
xlim([1, 2500])
ylim([-d d]);

plot(y_d, 'LineWidth', 2, ...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '\boldmath$y_d$')

for i = 1:N
    p(i) = plot(y1(i, :), 'color', colors2(i,:), ...
        'LineStyle', linestyle(i), ...
        'Marker', markers(i), ...
        'MarkerSize', 4, ...
        'LineWidth', 1.5, ...
        'DisplayName', ['\boldmath$y_{', num2str(i), '}$']);
end
yRange = ylim;
% xlabel('Time(k)', 'Interpreter', 'latex');
% ylabel('Output without Attacks', 'Interpreter', 'latex');
% lgd2=legend('Interpreter', 'latex', 'Location', positions2);
% lgd2.FontSize =6;
grid on
ax_sub1 = gca;
set(ax_sub1, 'Box', 'on'); % 添加边框



sub3 = subplot(2,2,4);
title('(d)')
% sub3.Position = subplot_positions4;

hold on
xlim([1, 2500])
ylim([-d d]);
for i = 1:T-1
    if DoS(i) == 0
        color = [0.95 0.95 0.95]; % 灰色
        patch([i, i+1, i+1, i], ...
            [yRange(1), yRange(1), yRange(2), yRange(2)], ...
            color, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

plot(y_d, 'LineWidth', 2, ...
    'Color', [0.8 0.8 0.8], ...
    'DisplayName', '\boldmath$y_d$')
for i = 1:N
    p(i) = plot(y(i, :), 'color', colors2(i,:), ...
        'LineStyle', linestyle(i), ...
        'Marker', markers(i), ...
        'MarkerSize', 4, ...
        'LineWidth', 1.5, ...
        'DisplayName', ['\boldmath$y_{', num2str(i), '}$']);
end
% xlabel('Time(k)', 'Interpreter', 'latex');
% ylabel('Output under Attacks', 'Interpreter', 'latex');
% lgd3=legend('Interpreter', 'latex', 'Location', positions3);
% lgd3.FontSize =6;
grid on
ax_sub2 = gca;
set(ax_sub2, 'Box', 'on'); % 添加边框

set(ax_sub2, 'Layer', 'top');

xlim([1, 2500])
subplot(2,2,2);
% sub4.Position = subplot_positions2;

hold on
for i = 1:N
    pr(i) = plot(r(i, :), 'color', colors2(i,:), ...
        'LineStyle', linestyle(i), ...
        'Marker', markers(i), ...
        'MarkerSize', 4, ...
        'LineWidth', 1.5, ...
        'DisplayName', ['\boldmath$agent_{', num2str(i), '}$']);
end
% xlabel('Time(k)', 'Interpreter', 'latex');
% ylabel('Event Triggering Function Value', 'Interpreter', 'latex');
lgd4 =legend('Interpreter', 'latex', 'Location', positions4);
lgd4.FontSize =6;
title('(b)')

grid on
ax_sub3 = gca;
set(ax_sub3, 'Box', 'on'); % 添加边框

xlim([1, 2500])





set(fig1, 'Units', 'inches');
pos = get(fig1, 'Position');
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'PaperSize', [pos(3), pos(4)]);
set(fig1, 'PaperPosition', [0 0 pos(3) pos(4)]);



set(fig2, 'Units', 'inches');
pos2 = get(fig2, 'Position');
set(fig2, 'PaperUnits', 'inches');
set(fig2, 'PaperSize', [pos2(3), pos2(4)]);
set(fig2, 'PaperPosition', [0 0 pos2(3) pos2(4)]);
%
% set(fig3, 'Units', 'inches');
% pos3 = get(fig3, 'Position');
% set(fig3, 'PaperUnits', 'inches');
% set(fig3, 'PaperSize', [pos3(3), pos3(4)]);
% set(fig3, 'PaperPosition', [0 0 pos3(3) pos3(4)]);
%
print(fig1, 'C:\Users\setup\OneDrive\Documents\Paper_4 2_15_2024\figure\Figure_2', '-dpdf')
print(fig2, 'C:\Users\setup\OneDrive\Documents\Paper_4 2_15_2024\figure\Figure_3', '-dpdf')
% print(fig3, 'C:\Users\setup\OneDrive\Documents\Paper_11_19_2023\code\fig6', '-dpdf')
% % W = diag([1,1,-1,1,1,1,-1,1]);
% % S = diag([3,3,4,3,3,3,4,3]);
% % A = [0,1,-1,0,0,0,0,0;
% %      0,0,0,1,0,0,0,0;
% %      0,0,0,0,0,0,0,0;
% %      0,0,0,0,1,0,0,0;
% %      0,0,0,0,0,1,0,0;
% %      0,0,0,0,0,0,0,1;
% %      0,0,0,0,0,0,0,0;
% %      0,0,0,0,0,0,-1,0];
% % absA = W*A*W;
% % D_absA = diag(absA*ones(8, 1));
% % L_absA = D_absA - absA;
% % D_hat = diag((S\A*S)*ones(8, 1));
% % L_absG = D_hat - S\A*S;
% % G = diag([1,0,0,0,0,1,0,0]);
% % varPsi2 = W*S*L_absG+G;
