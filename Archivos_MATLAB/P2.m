load datos_preg2.mat;
% Promedio Normal
normal_mean_erp = mean(EEGSegments_ERP,2);
plot(EEGTime_ERP, normal_mean_erp)
title('Promediación Normal de trials de ERP_ data')
xlabel('Tiempo (milisegundos)')
ylabel('Potencial evocado auditivo cortical (ERP_ data)')
figure()
normal_mean_vep = mean(EEGSegments_VEP,2);
plot(EEGTime_VEP, normal_mean_vep)
title('Promediación Normal de trials de SSVEP_ data')
xlabel('Tiempo (milisegundos)')
ylabel('Potencial evocado visual de estado estable (registro SSVEP_ data)')

x1_erp = 200;
x2_erp = 306;
x1_vep = 277;
x2_vep = 1563;
signal_erp = normal_mean_erp(x1_erp:x2_erp);
signal_vep = normal_mean_vep(x1_vep:x2_vep);
noise_erp = cat(1,normal_mean_erp(1:x1_erp - 1),normal_mean_erp(x2_erp+1:end));
noise_vep = cat(1,normal_mean_vep(1:x1_vep - 1),normal_mean_vep(x2_vep+1:end) );

snr_normal_erp = mean(signal_erp.^2)/mean(noise_erp.^2);
snr_normal_vep = mean(signal_vep.^2)/mean(noise_vep.^2);

%Promedio Ordenado
[~ ,ntrials_ERP] = size(EEGSegments_ERP);
[~ ,ntrials_VEP] = size(EEGSegments_VEP);

to_sort_ERP = zeros(1,ntrials_ERP);
for i = 1:ntrials_ERP
    to_sort_ERP(i) = var(EEGSegments_ERP(i,:));
end

[~,ind] = sort(to_sort_ERP);
sorted_ERP = EEGSegments_ERP(:,ind);
figure()
sorted_mean_erp = mean(sorted_ERP(:,1:round(ntrials_ERP/4)),2);
plot(EEGTime_ERP, sorted_mean_erp)
title('Promediación Ordenada de trials de ERP_ data')
xlabel('Tiempo (milisegundos)')
ylabel('Potencial evocado auditivo cortical (ERP_ data)')

signal_erp_sorted = sorted_mean_erp(x1_erp:x2_erp);
noise_erp_sorted = cat(1,sorted_mean_erp(1:x1_erp - 1),sorted_mean_erp(x2_erp+1:end));
snr_sorted_erp = mean(signal_erp_sorted.^2)/mean(noise_erp_sorted.^2);

% Promedio Ponderado
var_VEP = zeros(1,ntrials_VEP);
for i = 1:ntrials_VEP
    var_VEP(i) = var(EEGSegments_VEP(i,:));
end
sum_inv_var = sum(var_VEP.^-1).^-1;
pesos = var_VEP.^-1 * sum_inv_var;
pond_mean_vep = mean(EEGSegments_VEP.*repmat(pesos, [length(EEGTime_VEP) 1]), 2);
figure()
plot(EEGTime_VEP, pond_mean_vep)
title('Promediación Ponderada de trials de SSVEP_ data')
xlabel('Tiempo (milisegundos)')
ylabel('Potencial evocado visual de estado estable (registro SSVEP_ data)')

signal_vep_pond = pond_mean_vep(x1_vep:x2_vep);
noise_vep_pond = cat(1,pond_mean_vep(1:x1_vep - 1),pond_mean_vep(x2_vep+1:end));
snr_pond_vep = mean(signal_vep_pond.^2)/mean(noise_vep_pond.^2);

% Representacion Esprectral
sr_erp = 512;
sr_vep = 256;

L = length(EEGTime_ERP);
Y = fft(normal_mean_erp);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = sr_erp*(0:(L/2))/L;
figure()
plot(f, P1)
title('Espectrograma(FFT) de la promediación de ERP_ data')
xlabel('Frecuencia(Hz)')
ylabel('Potencial evocado auditivo cortical (ERP_ data)')

L = length(EEGTime_VEP);
Y = fft(normal_mean_vep);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = sr_erp*(0:(L/2))/L;
figure()
plot(f, P1)
title('Espectrograma(FFT) de la promediación de SSVEP_ data')
xlabel('Frecuencia(Hz)')
ylabel('Potencial evocado visual de estado estable (registro SSVEP_ data)')

figure()
cwt(double(normal_mean_erp),sr_erp)
figure()
cwt(double(normal_mean_vep),sr_vep)