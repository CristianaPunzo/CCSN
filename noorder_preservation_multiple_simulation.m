clear;
close all;
clc;

rng(6);  % Fissa il seme per riproducibilità

num_simulations = 100;  % Numero di simulazioni
T = 100;  % Numero di passi temporali
n = 5;  % Numero di agenti

coverage_cost_all = zeros(num_simulations, T);  % Matrice per salvare tutti i costi di coverage

% Valori nominali dei parametri (fissi per tutte le simulazioni)
lambda = [2.4, 2.8, 3, 4.2, 3.6]';
delta = [1.2, 4, 1.6, 3, 1.8]';
eta = [5, 6, 2, 4, 3]' * 0.001;
omega = [1.5, 12, 8, 0.5, 21]';
phi = [1/6, 1/3, 1/2, 1/4, 1/5]' * pi;

delta_bar = (delta + delta') / 2;

for sim = 1:num_simulations
    fprintf('Esecuzione simulazione %d/%d...\n', sim, num_simulations);
    
    % --- Generazione posizione iniziale casuale ---
    q = zeros(n,T);
    q(:,1) = sort(rand(n,1) * 2 * pi);  % Posizioni iniziali casuali ordinate
    u = zeros(n, T);

    % Matrici di distanza e errore iniziale
    d = zeros(n,n,T);
    e = zeros(n,n,T);

    for k = 1:T
        for i = 1:n
            for j = 1:n
                e(i, j, k) = (1+3*sin(omega(i)*k+10*(j-i).*phi(i))) * delta(i) / 4;
            end
        end
    end

    d(:,:,1) = angular_distance(q(:,1), q(:,1)') + e(:,:,1);
    d_bar = zeros(n,n,T);
    d_bar(:,:,1) = (d(:,:,1)-d(:,:,1)')/2;

    d_underhat = zeros(n,n,T);
    d_hat = zeros(n,n,T);
    
    d_underhat(:,:,1) = d_bar(:,:,1) - delta_bar;
    d_hat(:,:,1) = d_bar(:,:,1) + delta_bar;

    % --- Simulazione ---
    for i = 1:T-1
        d_tilde = (d_hat(:,:,i) + d_underhat(:,:,i)) / 2;

        lambda_shift_left = circshift(lambda, -1);
        lambda_shift_right = circshift(lambda, 1);

        d_tilde_forward = diag(circshift(d_tilde, -1, 2));
        d_tilde_backward = diag(circshift(d_tilde, 1, 2));

        u_tilde = eta .* ((lambda_shift_right + lambda) .* d_tilde_forward + ...
            (lambda + lambda_shift_left) .* d_tilde_backward);
        u(:,i) = lambda .* sat(u_tilde);

        q(:,i+1) = q(:,i) + u(:,i);

        d(:,:,i+1) = angular_distance(q(:,i+1), q(:,i+1)') + e(:,:,i+1);
        d_bar(:,:,i+1) = (d(:,:,i+1)-d(:,:,i+1)')/2;

        d_underhat(:,:,i+1) = max(d_bar(:,:,i+1) - delta_bar, d_underhat(:,:,i) + u(:,i)' - u(:,i));
        d_hat(:,:,i+1) = min(d_bar(:,:,i+1) + delta_bar, d_hat(:,:,i) + u(:,i)' - u(:,i));
    end

    % --- Calcolo della funzione di costo coverage per questa simulazione ---
    T_values_time = zeros(1, T);
    for t = 1:T
        T_values_time(t) = coverage_cost(q(:,t), lambda);
    end
    coverage_cost_all(sim, :) = T_values_time;
end

% --- Calcolo di media e deviazione standard ---
coverage_mean = mean(coverage_cost_all, 1);
coverage_std = std(coverage_cost_all, 0, 1);
T_star = pi / sum(lambda);  % Valore teorico

% --- Plot della media e deviazione standard ---
figure;
hold on;
x = 1:T;

% Banda di incertezza (± deviazione standard)
fill([x, fliplr(x)], [coverage_mean + coverage_std, fliplr(coverage_mean - coverage_std)], ...
    [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4); % Colore azzurro chiaro

% Plot della media
plot(x, coverage_mean, 'b-', 'LineWidth', 2, 'DisplayName', 'Media');

% Linea di riferimento teorica
yline(T_star, '--r', 'T^*', 'FontWeight', 'bold', 'DisplayName', 'T^*');

% Miglioramenti estetici
xlabel('Tempo (step)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Coverage cost T', 'FontSize', 14, 'FontWeight', 'bold');
title('Coverage cost medio su 100 simulazioni', 'FontSize', 16, 'FontWeight', 'bold');
legend('show', 'Location', 'best');
grid on;
hold off;

% Salva la figura
saveas(gcf, 'coverage_cost_mean_std.svg');

fprintf('Simulazioni completate! Grafico della coverage cost salvato.\n');

%% --- Funzioni di supporto ---
function d = angular_distance(x,y)
d_ccw = mod(y - x, 2*pi);
d_cw = d_ccw - 2*pi;

id_forw = circshift(eye(5),1,2);
id_backw = circshift(eye(5),-1,2);

d = id_forw .* d_ccw + id_backw .* d_cw;
end

function T = coverage_cost(q_i,lambda_i)
% Calcolo della funzione T
num_points = 1000;  % Risoluzione della circonferenza
q = linspace(0, 2*pi, num_points);  % Punti sulla circonferenza
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
