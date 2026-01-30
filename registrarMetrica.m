function registrarMetrica(metodo, listaEventos, nombreArchivo)
% registrarMetrica - Guarda las métricas detectadas en un archivo Excel.
%
% Entradas:
%   metodo       : String con el nombre del método (ej. 'Directo', 'MUSIC')
%   listaEventos : Array de estructuras retornadas por calcularMetricas
%   nombreArchivo: Path completo del archivo .xlsx

if isempty(listaEventos)
    fprintf('>>> No hay eventos para registrar en Excel.\n');
    return;
end

% Convertir struct array a tabla
T = struct2table(listaEventos);

% (Código de adición de columna Metodo eliminado)

% Reordenar columnas para legibilidad
% Orden deseado: Metodo, Time, Freq, Validacion, Diferencia, Delta, P_max, P_noise
% Orden deseado (Sin Metodo ni Validacion)
desiredOrder = {'Tiempo_Max_Seg', 'Freq_Hz', ...
    'Diferencia_Potencia_dBm', 'Diferencia_Simple_dBm', ...
    'Delta_Estabilidad_dBm', ...
    'P_max_dBm', 'P_ruido_dBm', 'P_ruido_Mean_dBm', 'P_promedio_dBm'};

% Verificar que existan (por si acaso cambió la struct)
existingCols = T.Properties.VariableNames;
validCols = intersect(desiredOrder, existingCols, 'stable');

T_final = T(:, validCols);

% Escribir Excel
try
    writetable(T_final, nombreArchivo);
    fprintf('>>> Reporte Excel guardado: %s\n', nombreArchivo);
catch ME
    fprintf('>>> Error guardando Excel: %s\n', ME.message);
end
end
