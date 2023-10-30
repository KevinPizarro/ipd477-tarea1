%% PREGUNTA3
ccc; load('RawData.mat');

%% P3.a) Promediacion Clasica
T = 1/sampling_rate;
t_total = T*EEGPoints;
pts_per_epoch = double(EEGPoints)/(double(t_total)/60.0); %327/60 = 5.45 min. => EEGPoints/5.45 =61440 pts per epoch
time_vector = 0:double(T):(double(t_total)-double(T));

data_out = zeros(chan,pts_per_epoch);
for i = 1:68
   data_out(i,:) = promediacion(EEG_data(i,:),pts_per_epoch,EEGPoints); 
end

figure(1)%GRAFICA DE LOS CANALES 62 AL 68
for i = 62:68
    subplot(3,3,i-61)
    plot(time_vector,EEG_data(i,:));
    xlabel('Time [ms]');ylabel('Voltage [uV]')
    title(strcat('Registro: ', locations(i).labels));
end 

%GRAFICAS DE COMPARACION PARA CANALES 62 Y 63. ORIGINAL VS PROM. CLASICA
for i = 62:63
    figure(i-60)
    plot(time_vector,EEG_data(i,:),'b')
    hold on 
    plot(time_vector(1,1:pts_per_epoch),data_out(i,:),'r')
    legend(["Original","Classic mean"]);
    xlim([0 pts_per_epoch*T])
    xlabel('Time [ms]');ylabel('Voltage [uV]')
    title(strcat('Registro: ', locations(i).labels));
end


%% P3.b) Independant Component Analisis (ICA)
components = 30;
ICA = jader(EEG_data(1:64,:), components);%calculamos ICA con jader
winv = pinv(ICA);%pseudo inversa de matriz ICA

figure(4) %GRAFICA DE COMPONENTES
for i = 1:30
    subplot(5,6,i);
    plot(ICA(i,:))
    title(strcat("Component ",string(i)))
end

disp("Calculando correlaciones cruzadas...")
vector = [];%guarda todos los valores de la correlacion cruzada
for i = 65:68
    for j = 1:30
        [C1,lag1] = xcorr(winv(j,:),EEG_data(i,:), 'none'); %correlacion cruzada entre canales externos y las componentes
        vector = [vector C1];
    end
end

vector = vector/max(abs(vector));%normalización entre -1 y 1
DelComponents = [];%matriz para guardar los componentes a eliminar

counter = 0; %recorrera todo el vector
for i = 65:68
    for j = 1:30
        for k = 1:length(C1)
            counter = counter + 1;
            if (vector(1,counter) >= 0.8)%Para correlaciones mayores al 80% se borrara dicho componente
                DelComponents(counter) = j; %recuperamos el componente al que pertence el dato
            end
        end
    end
end

DelComponents = nonzeros(DelComponents);%eliminamos los 0
DelComponents = unique(DelComponents,'rows');%eliminamos datos/canales redundantes
ICA_MOD = ICA; %MANTENGO INTACTO ICA
ICA_RAND = ICA; %para comparacion
ICA_MANUAL = ICA; %eliminacion manual

ICA_MOD(DelComponents,:)=[];%eliminacion de componentes encontrados
act=ICA_MOD*EEG_data(1:64,:); 
Reconstruction = pinv(ICA_MOD) * act; %reconstruccion de la señal


%% Automatico
figure(5)
subplot(2,1,1);
plot(time_vector, EEG_data(62,:))
title(strcat('Registro: ', locations(62).labels), " con artefactos.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
subplot(2,1,2);
plot(time_vector, Reconstruction(62,:))
title(strcat('Registro: ', locations(62).labels), " sin artefactos.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])

%% Pseudo-Aleatorio
rand_components = round(rand(3,1)*20+10);%componentes aleatorios entre los componentes 10 y 30
ICA_RAND(rand_components,:)=[];%eliminacion aleatoria
act_rand=ICA_RAND*EEG_data(1:64,:); 
Reconstruction_rand = pinv(ICA_RAND) * act_rand;%reconstruccion de la señal

figure(6)
subplot(2,1,1);
plot(time_vector, EEG_data(62,:))
title(strcat('Registro: ', locations(62).labels), " con artefactos.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
subplot(2,1,2);
plot(time_vector, Reconstruction_rand(62,:))
title(strcat('Registro: ', locations(62).labels), " componentes aleatorios.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])

%% Manual
ICA_MANUAL([1 2 3 4 5],:)=[];%eliminacion manual. primeros 5 por defecto
act_manual=ICA_MANUAL*EEG_data(1:64,:); 
Reconstruction_manual = pinv(ICA_MANUAL) * act_manual;%reconstruccion de la señal

figure(7)
subplot(2,1,1);
plot(time_vector, EEG_data(62,:))
title(strcat('Registro: ', locations(62).labels), " con artefactos.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
subplot(2,1,2);
plot(time_vector, Reconstruction_manual(62,:))
title(strcat('Registro: ', locations(62).labels), " eliminación manual de primeros 5 componentes.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])

%% Promediacion sin artefactos
data_out2 = zeros(chan, pts_per_epoch);

for i = 1:64
   data_out2(i,:) = promediacion(Reconstruction(i,:),pts_per_epoch,EEGPoints); 
end

figure(8)
plot(time_vector(1,1:pts_per_epoch),data_out(62,1:pts_per_epoch),'r')
hold on
plot(time_vector(1,1:pts_per_epoch),data_out2(62,:),'b')
xlim([0 pts_per_epoch*T])
xlabel('Time [ms]');ylabel('Voltage [uV]')
title(strcat('Registro: ', locations(62).labels, " promediados."));
legend(["Sin artefactos","Con artefactos"])

%% P3.c) ICA. Efecto de los componentes eliminados HARDCODED
ICA_MANUAL_1 = ICA; ICA_MANUAL_2 = ICA; ICA_MANUAL_3 = ICA; ICA_MANUAL_4 = ICA;
ICA_MANUAL_1([1 2 3 4 5],:)=[];%eliminacion manual. primeros 5
ICA_MANUAL_2([1 2 3 4 5 6 7 8 9 10],:)=[];%eliminacion manual. primeros 10
ICA_MANUAL_3([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15],:)=[];%eliminacion manual. primeros 15
ICA_MANUAL_4([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20],:)=[];%eliminacion manual. primeros 20
act_manual_1=ICA_MANUAL_1*EEG_data(1:64,:); 
act_manual_2=ICA_MANUAL_2*EEG_data(1:64,:); 
act_manual_3=ICA_MANUAL_3*EEG_data(1:64,:); 
act_manual_4=ICA_MANUAL_4*EEG_data(1:64,:); 
Reconstruction_manual_1 = pinv(ICA_MANUAL_1) * act_manual_1;%reconstruccion de la señal
Reconstruction_manual_2 = pinv(ICA_MANUAL_2) * act_manual_2;%reconstruccion de la señal
Reconstruction_manual_3 = pinv(ICA_MANUAL_3) * act_manual_3;%reconstruccion de la señal
Reconstruction_manual_4 = pinv(ICA_MANUAL_4) * act_manual_4;%reconstruccion de la señal

%GRAFICOS DE COMPARACION
figure(9)
subplot(2,2,1);%5 COMPONENTES
plot(time_vector, EEG_data(62,:),'b')
title(strcat('Registro: ', locations(62).labels), " con artefactos.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
hold on
plot(time_vector, Reconstruction_manual_1(62,:),'r')
title(strcat('Registro: ', locations(62).labels), " eliminación manual de primeros 5 componentes.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
legend(["Original", "25 componentes"])

subplot(2,2,2);%10 COMPONENTES
plot(time_vector, EEG_data(62,:),'b')
title(strcat('Registro: ', locations(62).labels), " con artefactos.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
hold on
plot(time_vector, Reconstruction_manual_2(62,:),'r')
title(strcat('Registro: ', locations(62).labels), " eliminación manual de primeros 10 componentes.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
legend(["Original", "20 componentes"])

subplot(2,2,3);%15 COMPONENTES
plot(time_vector, EEG_data(62,:),'b')
title(strcat('Registro: ', locations(62).labels), " con artefactos.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
hold on
plot(time_vector, Reconstruction_manual_3(62,:),'r')
title(strcat('Registro: ', locations(62).labels), " eliminación manual de primeros 15 componentes.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
legend(["Original", "15 componentes"])

subplot(2,2,4);%20 COMPONENTES
plot(time_vector, EEG_data(62,:),'b')
title(strcat('Registro: ', locations(62).labels), " con artefactos.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
hold on
plot(time_vector, Reconstruction_manual_4(62,:),'r')
title(strcat('Registro: ', locations(62).labels), " eliminación manual de primeros 20 componentes.")
xlabel('Time [ms]');ylabel('Voltage [uV]')
ylim([-280 200]); xlim([0 60])
legend(["Original", "10 componentes"])

%% FUNCIONES
function val_promedio = promediacion(datos,VENTANA,largo)
    total = 0;
    cont = 0;
    j = 1;
    
    while j+VENTANA-1<=largo
        cont = cont+1;
        
        
        ventana_i=datos(j:j+VENTANA-1);
        total = ventana_i+total;
        
        j=j+VENTANA;
    end
    val_promedio = total/cont;
end