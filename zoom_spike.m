clear; close all; clc;

% Carica i dati
data_folder = 'dati';
file_eta_wp = 'risultati_eta_wp.mat';

% Controllo se il file esiste
if exist(fullfile(data_folder, file_eta_wp), 'file')
    data_eta_wp = load(fullfile(data_folder, file_eta_wp));
else
    error('Il file %s non è stato trovato. Controlla la cartella dati.', file_eta_wp);
end

% Estrai i dati
T = size(data_eta_wp.u, 2); % Numero totale di passi temporali
u_wp_eta = data_eta_wp.u;    % Azioni di controllo
time = 1:T;                 % Vettore del tempo

% Intervallo di zoom
zoom_start = 1;
zoom_end = 20;

%% ---- Grafico dell'input di controllo con zoom ----
figure;
hold on;
plot(time, u_wp_eta', 'LineWidth', 1.5); % Plotta tutti gli agenti
xlabel('\textbf{Time (steps)}', 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$u_i$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('\textbf{Control Input}', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend(arrayfun(@(x) sprintf('$u_{%d}$', x), 1:size(u_wp_eta,1), 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location', 'best');
grid on;

% Troviamo i limiti per la regione di zoom
y_min = min(u_wp_eta(:, zoom_start:zoom_end), [], 'all');
y_max = max(u_wp_eta(:, zoom_start:zoom_end), [], 'all');

% Disegna un rettangolo nel grafico principale
rectangle('Position', [zoom_start, y_min, (zoom_end - zoom_start), (y_max - y_min)], ...
    'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '--');

% Creazione della finestra di zoom
axes_zoom = axes('Position', [0.65, 0.2, 0.25, 0.25]); % Posizione del riquadro di zoom
box on;
hold on;

% Plotta il controllo solo nella finestra di zoom
plot(time(zoom_start:zoom_end), u_wp_eta(:, zoom_start:zoom_end)', 'LineWidth', 1.5);
xlim([zoom_start, zoom_end]);
ylim([y_min - 0.01, y_max + 0.01]); % Adatta limiti asse Y
xlabel('Zoom Time', 'FontSize', 8);
ylabel('$u_i$', 'FontSize', 10, 'Interpreter', 'latex');
title('Zoomed Section', 'FontSize', 10);

hold off;

% Salva la figura
saveas(gcf, fullfile('immagini', 'control_input_wp_eta_zoom.svg'));
