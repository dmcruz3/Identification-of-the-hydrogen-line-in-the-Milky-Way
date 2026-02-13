function metrics = calcularMetricas(f, P_lin, t, config)
% calcularMetricas - Calcula métricas robustas para detección de HI.
%
% Entradas:
%   f      : Vector de frecuencias (Hz)
%   P_lin  : Matriz o Vector de Potencia LINEAL (mW o V^2).
%            NO dBm.
%            Si es matriz: [Frec x Tiempo]
%   t      : Vector de tiempo (s)
%
% Salida:
%   metrics : Struct con métricas en dBm (convertidas al final).

metrics = [];

% 0. Validaciones básicas
if isempty(P_lin) || isempty(f), return; end

% Asegurar dimensiones
if ismatrix(P_lin) && size(P_lin, 2) > 1
    % Perfil para detección de RUIDO: Promedio Lineal (Energía)
    P_perfil_avg_lin = mean(P_lin, 2);

    % Perfil para detección de PICO: Máximo Instantáneo
    P_perfil_max_lin = max(P_lin, [], 2);

    [~, idx_t_max] = max(P_lin, [], 2);
else
    P_perfil_avg_lin = P_lin(:);
    P_perfil_max_lin = P_lin(:);
    idx_t_max = ones(size(P_lin));
end

% 1. Encontrar Pico Máximo (Usando Perfil MAXIMO)
[P_max_val_lin, idx_peak] = max(P_perfil_max_lin);
freq_peak = f(idx_peak);

% Estimación robusta del ruido (Usando Perfil PROMEDIO para suavizar ruido)
ancho_excl = 20e3;
mask_noise = abs(f - freq_peak) > ancho_excl;

if any(mask_noise)
    P_ruido_val_lin = median(P_perfil_avg_lin(mask_noise));
else
    P_ruido_val_lin = median(P_perfil_avg_lin);
end

% Estimación Simple (Media Lineal)
P_ruido_mean_val_lin = mean(P_perfil_avg_lin);

% 2. Métrica Principal: RELACIÓN (SNR) Lineal
% Diferencia = P_max / P_ruido (en lineal es división)
SNR_lin = P_max_val_lin / P_ruido_val_lin;
SNR_Simple_lin = P_max_val_lin / P_ruido_mean_val_lin;

% 3. Análisis de Estabilidad
Delta_Estabilidad_dB = 0;
tiempo_evento = t(1);

if ismatrix(P_lin) && size(P_lin, 2) > 1
    traza_temporal_lin = P_lin(idx_peak, :);
    % Delta en dB: 10*log10(Max_lin / Min_lin)
    Delta_Estabilidad_dB = 10*log10(max(traza_temporal_lin) / min(traza_temporal_lin));
    tiempo_evento = t(idx_t_max(idx_peak));
end

% 4. Empaquetar Resultado (Convertir todo a dBm para reporte)
m_struct = struct();
m_struct.Freq_Hz = freq_peak;
m_struct.Tiempo_Max_Seg = tiempo_evento;

m_struct.P_max_dBm = 10*log10(P_max_val_lin);
m_struct.P_ruido_dBm = 10*log10(P_ruido_val_lin);
m_struct.P_ruido_Mean_dBm = 10*log10(P_ruido_mean_val_lin);
m_struct.P_promedio_dBm = 10*log10(mean(P_perfil_avg_lin));

% Diferencia en dB = 10*log10(SNR_lineal)
m_struct.Diferencia_Potencia_dBm = 10*log10(SNR_lin);
m_struct.Diferencia_Simple_dBm = 10*log10(SNR_Simple_lin);

m_struct.Delta_Estabilidad_dBm = Delta_Estabilidad_dB;

% -- Valores Lineales (NUEVOS - Pedido Ingeniero) --
m_struct.P_max_Lineal = P_max_val_lin;
m_struct.P_ruido_Lineal = P_ruido_val_lin;
m_struct.SNR_Lineal = SNR_lin; % Relación (veces)

metrics = [metrics, m_struct];
end
