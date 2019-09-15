load('Names.mat');
load('Dynamic.mat');
load('Static.mat');


% [Dynamic, Names] = Leer_archivo_walking(Pac);
% [Static, Names_standing]  = Leer_archivo_standing(Pac);
% 
% Medidas = Leer_archivo_medidas(Pac);
% 
% num_inicial= Dynamic(1,1);
% num_final= Dynamic(length(Dynamic(:,1)),1) + 1;

%Cargar datos al workspace antes de iniciar
%Establecer las medidas antropometricas del paciente, peso y talla


num_inicial=316;
num_final=709;



% num_inicial_rfp=524;
% num_final_rfp=623;
% num_inicial_lfp=445;
% num_final_lfp=542;


weight=64;

% height = Medidas(1);
% 
% R_mall_distance=Medidas(9);
% R_epy_distance=Medidas(7);
% R_pelvis_depth=Medidas(3);
% R_real_leg_length=Medidas(5);
% 
% L_mall_distance=Medidas(10);
% L_epy_distance=Medidas(8);
% L_pelvis_depth=Medidas(4);
% L_real_leg_length=Medidas(6);
% 
% Asis_distance=Medidas(2);


height=1.58;

R_mall_distance=0.07;
R_epy_distance=0.08;
R_pelvis_depth=0.155;
R_real_leg_length=0.805;

L_mall_distance=0.07;
L_epy_distance=0.08;
L_pelvis_depth=0.145;
L_real_leg_length=0.805;

Asis_distance=0.23;

Radius_marker=0.015;


[Frames, name_data]=size(Dynamic);

%Eventos

% Eventos = Leer_archivo_eventos(Pac);
% 
% r1_on=Eventos(1,2)*100+1-num_inicial;
% r2_on=Eventos(2,2)*100+1-num_inicial;
% r_off=Eventos(1,3)*100+1-num_inicial;
% l1_on=Eventos(1,4)*100+1-num_inicial;
% l2_on=Eventos(2,4)*100+1-num_inicial;
% l_off=Eventos(1,5)*100+1-num_inicial;

r1_on=521+2-num_inicial;
r2_on=668+2-num_inicial;
r_off=620+2-num_inicial;
l1_on=442+2-num_inicial;
l2_on=601+2-num_inicial;
l_off=544+2-num_inicial;


%Eje global LAB

Eje_Global=[1 0 0; 0 1 0; 0 0 1];

%Datos para el filtro

f_muestreo=100;
f_corte=6;
orden=2;

%Interpolación de los datos de las cámaras optoelectrónicas no debemos
%interpolar los datos de las plataformas y filtrado

% estan_plataformas=Find_name(Names,'r gr.X');
% 
% if isnan(estan_plataformas)
%     valor_final=name_data;
% else
%    valor_final=estan_plataformas-1; 
% end
% 
% 
% Dynamic_int=[];
% Dynamic_int_fil=[];
% %Dynamic_int(:,1)=Dynamic(num_inicial:num_final,1);
% %Dynamic_int(:,2)=Dynamic(num_inicial:num_final,2);
% %Dynamic_int_fil(:,1)=Dynamic(num_inicial:num_final,1);
% %Dynamic_int_fil(:,2)=Dynamic(num_inicial:num_final,2);
% 
% Dynamic_int(:,1)=Dynamic(:,1);
% Dynamic_int(:,2)=Dynamic(:,2);
% Dynamic_int_fil(:,1)=Dynamic(:,1);
% Dynamic_int_fil(:,2)=Dynamic(:,2);
% 
% 
% for x=1:length(Dynamic(1,:))
%     
%     %Dynamic_int(:,x)=Interpolacion_cubica(Dynamic(num_inicial:num_final,1),Dynamic(num_inicial:num_final,x),0);
%     Dynamic_int(:,x)=Interpolacion_cubica(Dynamic(:,1),Dynamic(:,x),0);
%     Dynamic_int_fil(:,x)=Butterworth_filter(Dynamic_int(:,x),f_muestreo, f_corte, orden);
%     
% end




%% Standing

%Corrección de los asis por medio de las medidas antropometricas
[SR_Asis2,SL_Asis2]=Asis_correction(Names,Static,Asis_distance);