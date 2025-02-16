clear
close all
clc

num_simulations = 100;
T = 100;
n = 5;

lambda = [2.4, 2.8, 3, 4.2, 3.6]';
delta = [0.32, 0.24, 0.45, 0.18, 0.36]';
eta = [0.010, 0.007, 0.008, 0.006, 0.009]';
omega = [1.5, 12, 8, 0.5, 21]';
phi=[1/6, 1/3, 1/2, 1/4, 1/5]'*pi;

T_values_all = zeros(num_simulations, T);

for sim = 1:num_simulations
    rng(sim)
    q = zeros(n,T);
    q(:,1) = sort(rand(n,1)*2*pi);
    u = zeros(n,T);
    d_underhat = zeros(n,n,T);
    d_hat = zeros(n,n,T);
    delta_bar = zeros(n, n);
    
    for i = 1:n
        for j = 1:n
            delta_bar(i, j) = (delta(i) + delta(j)) / 2;
        end
    end

    e = zeros(n,n,T);
    for k = 1:T
        for i=1:n
            for j=1:n
                e(i, j, k) = (1+3*sin(omega(i)*k+10*(j-i).*phi(i)))*delta(i)/4;
            end
        end
    end

    d = zeros(n,n,T);
    d(:,:,1) = angular_distance(q(:,1),q(:,1)') + e(:,:,1);
    d_bar = zeros(n,n,T);
    d_bar(:,:,1) = (d(:,:,1)-d(:,:,1)')/2;
    d_underhat(:,:,1) = d_bar(:,:,1) - delta_bar;
    d_hat(:,:,1) = d_bar(:,:,1) + delta_bar;

    for i = 1:T-1
        d_tilde = (d_hat(:,:,i) + d_underhat(:,:,i)) / 2;
        lambda_shift_left = circshift(lambda, -1);  
        lambda_shift_right = circshift(lambda, 1);  
        d_underhat_i = d_underhat(:,:,i);
        d_underhat_forward = diag(circshift(d_underhat_i, -1, 2)); 
        d_underhat_backward = diag(circshift(d_underhat_i, 1, 2)); 

        u_tilde = eta .* ((lambda_shift_right + lambda) .* d_underhat_forward + ...
            (lambda + lambda_shift_left) .* d_underhat_backward);
        u(:,i) = lambda .* sat(sigma(u_tilde));
        q(:,i+1) = q(:,i) + u(:,i);
        d(:,:,i+1) = angular_distance(q(:,i+1),q(:,i+1)') + e(:,:,i+1);
        d_bar(:,:,i+1) = (d(:,:,i+1)-d(:,:,i+1)')/2;
        d_underhat(:,:,i+1) = max(d_bar(:,:,i+1) - delta_bar, d_underhat(:,:,i) + u(:,i)' - u(:,i));
        d_hat(:,:,i+1) = min(d_bar(:,:,i+1) + delta_bar, d_hat(:,:,i) + u(:,i)' - u(:,i));
    end

    for t = 1:T
        T_values_all(sim, t) = coverage_cost(q(:,t), lambda);
    end
end

T_star = pi / sum(lambda);
T_mean = mean(T_values_all, 1);
T_std = std(T_values_all, 0, 1);

figure;
hold on;
x = 1:T;
fill([x, fliplr(x)], [T_mean + T_std, fliplr(T_mean - T_std)], ...
    [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(x, T_mean, 'b-', 'LineWidth', 2);
yline(T_star, '--r', 'T^*', 'FontWeight', 'bold');
xlabel('\textbf{Time (steps)}', 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$T$',  'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('\textbf{Coverage Cost Function over 100 simulations}', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend({'Mean', 'T'}, 'Location', 'best');
grid on;
hold off;
saveas(gcf, 'cost_function_wp_mean_std_wp.svg');

% Determina gli ultimi valori da considerare come regime (ultimo 10% dei passi temporali)
final_portion = round(0.1 * T); 
steady_state_mean = mean(T_mean(end-final_portion+1:end));
steady_state_std = mean(T_std(end-final_portion+1:end));

% Stampa i risultati
fprintf('Valore a regime della media: %.6f\n', steady_state_mean);
fprintf('Valore a regime della deviazione standard: %.6f\n', steady_state_std);



function d = angular_distance(x,y)
d_ccw = mod(y - x, 2*pi);
d_cw = d_ccw - 2*pi;
id_forw = circshift(eye(5),1,2);
id_backw = circshift(eye(5),-1,2);
d = id_forw .* d_ccw + id_backw .* d_cw;
end

function T = coverage_cost(q_i,lambda_i)
num_points = 1000;
q = linspace(0, 2*pi, num_points);
T_values = zeros(1, num_points);
for k = 1:num_points
    dist = mod(q(k) - q_i, 2*pi);
    dist = min(dist, 2*pi - dist);
    min_vals = dist ./ lambda_i;
    T_values(k) = min(min_vals);
end
T = max(T_values);
end

function y = sat(u)
y = sign(u).*min(1,abs(u));
end

function s = sigma(u)
s = max(u,0);
end
