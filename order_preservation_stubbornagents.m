clear
close all
clc

rng(6)

n = 6;
T = 200;

q = zeros(n,T);
q(:,1) = sort(rand(n,1)*2*pi);
u = zeros(n,T);
lambda = [0.01, 2.4, 2.8, 0.01, 4.2, 3.6]';
delta = [0.1, 0.32, 0.24, 0.45, 0.18, 0.36]';
d_underhat = zeros(n,n,T);
d_hat = zeros(n,n,T);
delta_bar = zeros(n, n);
eta = [0, 0.010, 0.007, 0.008, 0.006, 0.009]'; %valori nominali
omega = [5, 1.5, 12, 8, 0.5, 21]';
phi = [1/2, 1/6, 1/3, 1/2, 1/4, 1/5]'*pi;

for i = 1:n
    for j = 1:n
        delta_bar(i, j) = (delta(i) + delta(j)) / 2;
    end
end

e = zeros(n,n,T);
for k = 1:T
    for i = 1:n
        for j = 1:n
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

% Genera una palette di colori distinti per ogni agente
colors = lines(n);  

% Visualizzazione posizione iniziale
figure;
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'k', 'LineWidth', 1.5);
hold on;
x_init = cos(q(:,1));
y_init = sin(q(:,1));
for j = 1:n
    scatter(x_init(j), y_init(j), 100, colors(j,:), 'filled');
    text(1.1*x_init(j), 1.1*y_init(j), num2str(j), 'FontSize', 12, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Color', colors(j,:), 'Interpreter', 'latex');
end
axis equal;
xlim([-1.2 1.2]); ylim([-1.2 1.2]);
title('\textbf{Initial Position}', 'Interpreter', 'latex');
hold off;
saveas(gcf, fullfile('immagini', 'posizione_iniziale_stubborn_wp.svg'));

% Inizializza VideoWriter per salvare l'animazione in formato MP4
video_filename = fullfile('immagini', 'stubborn_wp.mp4');
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 10;  % Regola il frame rate secondo le necessità
open(v);

figure;
for i = 1:T-1
    d_tilde = (d_hat(:,:,i) + d_underhat(:,:,i)) / 2;
    
    lambda_shift_left = circshift(lambda, -1);  
    lambda_shift_right = circshift(lambda, 1);  
    
    d_underhat_i = d_underhat(:,:,i);
    % Calcola d_underhat_forward e d_underhat_backward
    d_underhat_forward = diag(circshift(d_underhat_i, -1, 2)); 
    d_underhat_backward = diag(circshift(d_underhat_i, 1, 2)); 
    
    u_tilde = eta .* ((lambda_shift_right + lambda) .* d_underhat_forward + ...
        (lambda + lambda_shift_left) .* d_underhat_backward);
    u(:,i) = lambda .* sat(sigma(u_tilde));
    
    % Imposta gli agenti "stubborn" (in questo caso 1 e 4)
    u(1,i) = lambda(1);
    u(4,i) = lambda(4);
    
    q(:,i+1) = q(:,i) + u(:,i);
    
    % --- Visualizzazione aggiornata ---
    clf;
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k', 'LineWidth', 1.5); % Disegna la circonferenza
    hold on;
    
    % Converte le coordinate polari in cartesiane
    x = cos(q(:,i));
    y = sin(q(:,i));
    
    for j = 1:n
        scatter(x(j), y(j), 100, colors(j,:), 'filled');
    end
    
    % Etichettatura degli agenti
    x_text = 1.1 * x;
    y_text = 1.1 * y;
    for j = 1:n
        text(x_text(j), y_text(j), num2str(j), 'FontSize', 12, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Color', colors(j,:));
    end
    hold off;
    axis equal; 
    xlim([-1.2 1.2]); ylim([-1.2 1.2]); 
    title(sprintf('Step %d', i));
    drawnow;
    
    % Cattura il frame corrente e lo scrive nel video
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    d(:,:,i+1) = angular_distance(q(:,i+1),q(:,i+1)') + e(:,:,i+1);
    d_bar(:,:,i+1) = (d(:,:,i+1)-d(:,:,i+1)')/2;
    
    d_underhat(:,:,i+1) = max(d_bar(:,:,i+1) - delta_bar, d_underhat(:,:,i) + u(:,i)' - u(:,i));
    d_hat(:,:,i+1) = min(d_bar(:,:,i+1) + delta_bar, d_hat(:,:,i) + u(:,i)' - u(:,i));
end

% Chiude il VideoWriter
close(v);

% %% ---- Calcolo della funzione coverage T per tutti gli istanti di tempo ----
% T_values_time = zeros(1, T);
% for t = 1:T
%     T_values_time(t) = coverage_cost(q(:,t), lambda);
% end
% 
% T_star = pi / sum(lambda);
% 
% figure;
% plot(1:T, T_values_time - T_star, 'b-', 'LineWidth', 1.5);
% yline(0, '--k', 'T^*', 'FontWeight', 'bold');
% xlabel('\textbf{Time (steps)}', 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('$T - T^*$',  'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% title('\textbf{Coverage Cost Function}', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% grid on;
% grid minor;
% ylim([min(T_values_time - T_star) - 0.01, max(T_values_time - T_star) + 0.01]);
% xlim([0, T]);
% saveas(gcf, fullfile('immagini', 'cost_function_stubborn_wp.svg'));
% 
% %% ---- Posizione finale ----
% figure;
% plot(cos(theta), sin(theta), 'k', 'LineWidth', 1.5); % Disegna la circonferenza
% hold on;
% x_final = cos(q(:,end));
% y_final = sin(q(:,end));
% for j = 1:n
%     scatter(x_final(j), y_final(j), 100, colors(j,:), 'filled');
%     text(1.1*x_final(j), 1.1*y_final(j), num2str(j), 'FontSize', 12, ...
%          'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', colors(j,:));
% end
% axis equal;
% xlim([-1.2 1.2]); ylim([-1.2 1.2]);
% title('\textbf{Final Position}', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% hold off;
% saveas(gcf, fullfile('immagini', 'posizione_finale_stubborn_wp.svg'));
% 
% %% ---- Grafico delle posizioni nel tempo ----
% figure;
% hold on;
% for j = 1:n
%     plot(1:T, q(j,:), 'Color', colors(j,:), 'LineWidth', 1.5);
%     text(T, q(j,end), sprintf('q_{%d}', j), 'FontSize', 12, 'FontWeight', 'bold', ...
%          'Color', colors(j,:), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
% end
% xlabel('\textbf{Time (steps)}', 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('$q_i$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% title('\textbf{Position in time}', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% grid on;
% hold off;
% saveas(gcf, fullfile('immagini', 'position_stubborn_wp.svg'));
% 
% %% ---- Grafico del controllo nel tempo ----
% figure;
% plot(1:T, u', 'LineWidth', 1.5);
% xlabel('\textbf{Time (steps)}', 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'latex');
% ylabel('$u_i$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
% title('\textbf{Control Input}', 'FontSize', 14, 'Interpreter', 'latex');
% grid on;
% legend(arrayfun(@(x) sprintf('$u_{%d}$', x), 1:n, 'UniformOutput', false), 'Interpreter', 'latex');
% hold off;
% saveas(gcf, fullfile('immagini', 'control_input_stubborn_wp.svg'));
% 
% %% ---- Salvataggio dati ----
% if ~exist('dati', 'dir')
%     mkdir('dati');
% end
% save(fullfile('dati', 'risultati_stubborn_wp.mat'), 'T_values_time', 'T_star', 'u');

%% ---- Funzioni ausiliarie ----
function d = angular_distance(x,y)
    d_ccw = mod(y - x, 2*pi);
    d_cw = d_ccw - 2*pi;
    id_forw = circshift(eye(6),1,2);
    id_backw = circshift(eye(6),-1,2);
    d = id_forw .* d_ccw + id_backw .* d_cw;
end

function T = coverage_cost(q_i,lambda_i)
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

function s = sigma(u)
    s = max(u,0);
end
