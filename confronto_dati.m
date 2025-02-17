clear; close all; clc;

% Cartella contenente i file di dati
data_folder = 'dati';

% Lista dei file da confrontare (modifica con i tuoi file effettivi)
file_list = {'risultati_nominali.mat', 'risultati_pinning_wop.mat'}; % Assicurati che questi file esistano

% Inizializza variabili per la memorizzazione dei dati
T_values_all = {};
T_star_all = [];
u_all = {};
labels = {'Without Pinning', 'With Pinning'}; % Etichette della legenda

% Carica i dati dai file
for i = 1:length(file_list)
    file_path = fullfile(data_folder, file_list{i});
    if exist(file_path, 'file')
        data = load(file_path);
        T_values_all{i} = data.T_values_time;
        T_star_all(i) = data.T_star;
        u_all{i} = data.u;
    else
        error('File %s non trovato. Assicurati che sia nella cartella corretta.', file_list{i});
    end
end

% Numero di passi temporali
T = length(T_values_all{1});

%% ---- Confronto della funzione di costo T - T* ----
figure;
hold on;
plot(1:T, T_values_all{1} - T_star_all(1), 'b-', 'LineWidth', 1.5);
plot(1:T, T_values_all{2} - T_star_all(2), 'LineWidth', 1.5);
hold off;
xlabel('\textbf{Time (steps)}', 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$T - T^*$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('\textbf{Comparison of Coverage Cost Function}', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend(labels, 'Interpreter', 'latex', 'Location', 'Best');
grid on;
saveas(gcf, fullfile('immagini', 'comparison_cost_function_pinning_wop.svg'));

%% ---- Confronto delle azioni di controllo ----
figure;
hold on;
plot(1:T, mean(u_all{1}, 1), 'b-', 'LineWidth', 1.5); % Nominal
plot(1:T, mean(u_all{2}, 1), 'LineWidth', 1.5); % Non-Nominal
hold off;
xlabel('\textbf{Time (steps)}', 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$\bar{u}$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('\textbf{Comparison of Control Input}', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
legend(labels, 'Interpreter', 'latex', 'Location', 'Best');
grid on;
saveas(gcf, fullfile('immagini', 'comparison_control_input_pinning_wop.svg'));
