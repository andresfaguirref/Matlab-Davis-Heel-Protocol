function processed_signal=Butterworth_filter(datos, f_muestreo, f_corte, orden)

[num,den]=butter(orden,2*f_corte/f_muestreo, 'low');

processed_signal=filtfilt(num,den,datos);

end