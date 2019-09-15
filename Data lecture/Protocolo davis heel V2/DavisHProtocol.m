clear
clc

FData = {};
FNames = {};

for Pac = 1:7

    Pac

%% lectura de archivos

%Pnum = num2str(Pac);
%Nombre = ['Marcadores_Standing_P' Pnum '.txt'];


%%
       

load('Names.mat');
load('Dynamic.mat');
load('Static.mat');

Dynamic = [];
Names = {};
[Dynamic, Names] = Leer_archivo_walking(Pac);
[Static, Names_standing]  = Leer_archivo_standing(Pac);

Medidas = Leer_archivo_medidas(Pac);

num_inicial= Dynamic(1,1);
num_final= Dynamic(length(Dynamic(:,1)),1) + 1;

% Cargar datos al workspace antes de iniciar
% Establecer las medidas antropometricas del paciente, peso y talla


% num_inicial=316;
% num_final=709;

% num_inicial_rfp=524;
% num_final_rfp=623;
% num_inicial_lfp=445;
% num_final_lfp=542;


weight=64;

height = Medidas(1);

R_mall_distance=Medidas(9);
R_epy_distance=Medidas(7);
R_pelvis_depth=Medidas(3);
R_real_leg_length=Medidas(5);

L_mall_distance=Medidas(10);
L_epy_distance=Medidas(8);
L_pelvis_depth=Medidas(4);
L_real_leg_length=Medidas(6);

Asis_distance=Medidas(2);


% height=1.58;
% 
% R_mall_distance=0.07;
% R_epy_distance=0.08;
% R_pelvis_depth=0.155;
% R_real_leg_length=0.805;
% 
% L_mall_distance=0.07;
% L_epy_distance=0.08;
% L_pelvis_depth=0.145;
% L_real_leg_length=0.805;
% 
% Asis_distance=0.23;

Radius_marker=0.015;


[Frames, name_data]=size(Dynamic);

%Eventos

Eventos = Leer_archivo_eventos(Pac);

r1_on=int64(Eventos(1,1)*100+1-num_inicial);
r2_on=int64(Eventos(2,1)*100+1-num_inicial);
r_off=int64(Eventos(1,2)*100+1-num_inicial);
l1_on=int64(Eventos(1,3)*100+1-num_inicial);
l2_on=int64(Eventos(2,3)*100+1-num_inicial);
l_off=int64(Eventos(1,4)*100+1-num_inicial);

% r1_on=521+2-num_inicial;
% r2_on=668+2-num_inicial;
% r_off=620+2-num_inicial;
% l1_on=442+2-num_inicial;
% l2_on=601+2-num_inicial;
% l_off=544+2-num_inicial;


%Eje global LAB

Eje_Global=[1 0 0; 0 1 0; 0 0 1];

%Datos para el filtro

f_muestreo=100;
f_corte=6;
orden=2;

%Interpolación de los datos de las cámaras optoelectrónicas no debemos
%interpolar los datos de las plataformas y filtrado

estan_plataformas=Find_name(Names,'r gr.X');

if isnan(estan_plataformas)
    valor_final=name_data;
else
   valor_final=estan_plataformas-1; 
end


Dynamic_int=[];
Dynamic_int_fil=[];
%Dynamic_int(:,1)=Dynamic(num_inicial:num_final,1);
%Dynamic_int(:,2)=Dynamic(num_inicial:num_final,2);
%Dynamic_int_fil(:,1)=Dynamic(num_inicial:num_final,1);
%Dynamic_int_fil(:,2)=Dynamic(num_inicial:num_final,2);

Dynamic_int(:,1)=Dynamic(:,1);
Dynamic_int(:,2)=Dynamic(:,2);
Dynamic_int_fil(:,1)=Dynamic(:,1);
Dynamic_int_fil(:,2)=Dynamic(:,2);


for x=3:valor_final
    
    %Dynamic_int(:,x)=Interpolacion_cubica(Dynamic(num_inicial:num_final,1),Dynamic(num_inicial:num_final,x),0);
    Dynamic_int(:,x)=Interpolacion_cubica(Dynamic(:,1),Dynamic(:,x),0);
    Dynamic_int_fil(:,x)=Butterworth_filter(Dynamic_int(:,x),f_muestreo, f_corte, orden);
    
end









%% Standing

%Corrección de los asis por medio de las medidas antropometricas
[SR_Asis2,SL_Asis2]=Asis_correction(Names,Static,Asis_distance);

%Sistemas de referencia global para la marcha

RefS_global_marcha=Eje_Global_marcha(Names, Static);

%Sistemas de referencia técnicos y centros de rotación

%Rodilla

SRC_right_knee=SRight_Knee_RC(Names, Static, R_epy_distance);
SRC_left_knee=SLeft_Knee_RC(Names, Static, L_epy_distance);

%Tobillo

SRC_right_ankle=SRight_Ankle_RC(Names, Static, R_mall_distance);
SRC_left_ankle=SLeft_Ankle_RC(Names, Static, L_mall_distance);

%Metatarso

SRC_right_metatarsal=SRight_Metatarsal_RC(Names, Static, R_mall_distance);
SRC_left_metatarsal=SLeft_Metatarsal_RC(Names, Static, L_mall_distance);

%Pelvis

%Por metodo de Bell

SBell_RC_right_hip=SBell_RPelvis_RC(Names,SR_Asis2, SL_Asis2,Static);
SBell_RC_left_hip=SBell_LPelvis_RC(Names,SR_Asis2, SL_Asis2,Static);

%Con Paper Davis 1991

Leg_length=0.5*(L_real_leg_length+R_real_leg_length);
SDavis_RC_right_hip=SDavisH_RPelvis_RC(Names, Static, SR_Asis2, SL_Asis2, Leg_length);
SDavis_RC_left_hip=SDavisH_LPelvis_RC(Names, Static, SR_Asis2, SL_Asis2, Leg_length);



%Método del 25%

SPC25_RC_right_hip=SPerC25_RPelvis_RC(Names, Static);
SPC25_RC_left_hip=SPerC25_LPelvis_RC(Names, Static);

%Por ultimo método de BTS

RHPTP=RPelvis_Valor(Asis_distance, R_pelvis_depth, R_real_leg_length);
LHPTP=LPelvis_Valor(Asis_distance, L_pelvis_depth, L_real_leg_length);

[SBTS_RC_right_hip,SEje_Pelvis]=SBTS_RPelvis_RC(Names, SR_Asis2, SL_Asis2, Static, RHPTP);
SBTS_RC_left_hip=SBTS_LPelvis_RC(Names, SR_Asis2, SL_Asis2, Static, LHPTP);


%Sistema de referencia de el tronco 

[SRef_sys_Trunk,STrunk_RC]=RefS_Trunk(Names, Static);



%% Reference system foot progression

for x=1:length(SRC_right_ankle)
   
    [SRFP_RefS{x},SLFP_RefS{x}] =FP_RefS(RefS_global_marcha, SRC_right_ankle(x,:), SRC_right_metatarsal(x,:), SRC_left_ankle(x,:), SRC_left_metatarsal(x,:));
    
end

 


%% Ejes de referencia anatómicos

%Eje anatómico del Muslo por método de Bell

SRef_sys_RMuslo_Bell=SRefS_RMuslo(Names, Static, SBell_RC_right_hip, SRC_right_knee);

for x=1:length(SRef_sys_RMuslo_Bell)
    X_new = SRef_sys_RMuslo_Bell{x}(1,:);
    Y_new = SRef_sys_RMuslo_Bell{x}(3,:);
    Z_new = SRef_sys_RMuslo_Bell{x}(2,:).*-1;
    SRef_sys_RMuslo_Bell{x}=[X_new;Y_new;Z_new];    
end

SRef_sys_LMuslo_Bell=RefS_LMuslo(Names, Static, SBell_RC_left_hip, SRC_left_knee);

for x=1:length(SRef_sys_LMuslo_Bell)
    X_new = SRef_sys_LMuslo_Bell{x}(1,:);
    Y_new = SRef_sys_LMuslo_Bell{x}(3,:).*-1;
    Z_new = SRef_sys_LMuslo_Bell{x}(2,:);
    SRef_sys_LMuslo_Bell{x}=[X_new;Y_new;Z_new];    
end

%Eje anatómico del Muslo por método del 25%

SRef_sys_RMuslo_25=SRefS_RMuslo(Names, Static, SPC25_RC_right_hip, SRC_right_knee);

for x=1:length(SRef_sys_RMuslo_25)
    X_new = SRef_sys_RMuslo_25{x}(1,:);
    Y_new = SRef_sys_RMuslo_25{x}(3,:);
    Z_new = SRef_sys_RMuslo_25{x}(2,:).*-1;
    SRef_sys_RMuslo_25{x}=[X_new;Y_new;Z_new];    
end

SRef_sys_LMuslo_25=SRefS_LMuslo(Names, Static, SPC25_RC_left_hip, SRC_left_knee);

for x=1:length(SRef_sys_LMuslo_25)
    X_new = SRef_sys_LMuslo_25{x}(1,:);
    Y_new = SRef_sys_LMuslo_25{x}(3,:).*-1;
    Z_new = SRef_sys_LMuslo_25{x}(2,:);
    SRef_sys_LMuslo_25{x}=[X_new;Y_new;Z_new];    
end

%Eje anatómico del Muslo por método de Davis 1991

SRef_sys_RMuslo_Davis=SRefS_RMuslo(Names, Static, SDavis_RC_right_hip, SRC_right_knee);

for x=1:length(SRef_sys_RMuslo_Davis)
    X_new = SRef_sys_RMuslo_Davis{x}(1,:);
    Y_new = SRef_sys_RMuslo_Davis{x}(3,:);
    Z_new = SRef_sys_RMuslo_Davis{x}(2,:).*-1;
    SRef_sys_RMuslo_Davis{x}=[X_new;Y_new;Z_new];    
end

SRef_sys_LMuslo_Davis=SRefS_LMuslo(Names, Static, SDavis_RC_left_hip, SRC_left_knee);

for x=1:length(SRef_sys_LMuslo_Davis)
    X_new = SRef_sys_LMuslo_Davis{x}(1,:);
    Y_new = SRef_sys_LMuslo_Davis{x}(3,:).*-1;
    Z_new = SRef_sys_LMuslo_Davis{x}(2,:);
    SRef_sys_LMuslo_Davis{x}=[X_new;Y_new;Z_new];    
end


%Eje anatómico del Muslo por método de BTS

SRef_sys_RMuslo_BTS=SRefS_RMuslo(Names, Static, SBTS_RC_right_hip, SRC_right_knee);

for x=1:length(SRef_sys_RMuslo_BTS)
    X_new = SRef_sys_RMuslo_BTS{x}(1,:);
    Y_new = SRef_sys_RMuslo_BTS{x}(3,:);
    Z_new = SRef_sys_RMuslo_BTS{x}(2,:).*-1;
    SRef_sys_RMuslo_BTS{x}=[X_new;Y_new;Z_new];    
end

% for x=1:length(SRef_sys_RMuslo_BTS)
%    
%     X(x,:)=SRef_sys_RMuslo_BTS{x}(1,:);    
%        
% end
% 
% figure();
% subplot(3,1,1);
% plot(X(:,1));
% subplot(3,1,2);
% plot(X(:,2));
% subplot(3,1,3);
% plot(X(:,3));

SRef_sys_LMuslo_BTS=SRefS_LMuslo(Names, Static, SBTS_RC_left_hip, SRC_left_knee);

for x=1:length(SRef_sys_LMuslo_BTS)
    X_new = SRef_sys_LMuslo_BTS{x}(1,:);
    Y_new = SRef_sys_LMuslo_BTS{x}(3,:).*-1;
    Z_new = SRef_sys_LMuslo_BTS{x}(2,:);
    SRef_sys_LMuslo_BTS{x}=[X_new;Y_new;Z_new];    
end

%Eje anatómico de la Pierna

SRef_sys_RPierna=SRefS_RPierna(Names, Static, SRC_right_knee, SRC_right_ankle);
SRef_sys_LPierna=SRefS_LPierna(Names, Static, SRC_left_knee, SRC_left_ankle);

%Eje anatómico del pie 

SRef_sys_RPie_BTS=SRefS_RPie_BTS(Names, Static, SRC_right_ankle, SRC_right_metatarsal);
SRef_sys_LPie_BTS=SRefS_LPie_BTS(Names, Static, SRC_left_ankle, SRC_left_metatarsal);

%Eje de referencia de la pelvis

SRef_sys_Pelvis = {};
for x=1:length(SEje_Pelvis)
    X_new = SEje_Pelvis{x}(3,:);
    Y_new = SEje_Pelvis{x}(1,:);
    Z_new = SEje_Pelvis{x}(2,:);
    SRef_sys_Pelvis{x}=[X_new;Y_new;Z_new];    
end

%Eje de referencia de el tronco

for x=1:length(SRef_sys_Trunk)
    X_new = SRef_sys_Trunk{x}(3,:);
    Y_new = SRef_sys_Trunk{x}(2,:).*-1;
    Z_new = SRef_sys_Trunk{x}(1,:);
    SRef_sys_Trunk{x}=[X_new;Y_new;Z_new];    
end






%%
%Encontrar los centros de masa de los segmentos

%Foot segment

SRFoot_CoM=CoM(SRC_right_ankle,SRC_right_metatarsal,0.44);
SLFoot_CoM=CoM(SRC_left_ankle,SRC_left_metatarsal,0.44);

%Shank segment

SRShank_CoM=CoM(SRC_right_knee,SRC_right_ankle,0.42);
SLShank_CoM=CoM(SRC_left_knee,SRC_left_ankle,0.42);

%Thigh Segment

%With Bell method

SRThigh_CoM=CoM(SBell_RC_right_hip, SRC_right_knee,0.39);
SLThigh_CoM=CoM(SBell_RC_left_hip, SRC_left_knee,0.39);

%With 25% method

SR2Thigh_CoM=CoM(SPC25_RC_right_hip, SRC_right_knee,0.39);
SL2Thigh_CoM=CoM(SPC25_RC_left_hip, SRC_left_knee,0.39);

%with helen hayes (davis heel)

SR3Thigh_CoM=CoM(SDavis_RC_right_hip, SRC_right_knee,0.39);
SL3Thigh_CoM=CoM(SDavis_RC_left_hip, SRC_left_knee,0.39);

%With BTS method

SR4Thigh_CoM=CoM(SBTS_RC_right_hip, SRC_right_knee,0.39);
SL4Thigh_CoM=CoM(SBTS_RC_left_hip, SRC_left_knee,0.39);






%%

%Calcular ángulos del pie (dorsi y plantar flexión, inversión-eversión,
%rotación interna y externa)

%Con el método de BTS

for x=1:length(SRef_sys_RPierna)
    
SR_ankle_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_RPierna{x}, SRef_sys_RPie_BTS{x});
SR_ankle_angle(x,:)=SR_ankle_angle(x,:)+[0 0 90];
SR_ankle_angle(x,:)=SR_ankle_angle(x,:).*[-1, 1, -1];

end

SR_ankle_IntExt=mean(SR_ankle_angle(:,1));
SR_ankle_IE=mean(SR_ankle_angle(:,2));
SR_ankle_FE=mean(SR_ankle_angle(:,3));


for x=1:length(SRef_sys_LPierna)
    
SL_ankle_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_LPierna{x}, SRef_sys_LPie_BTS{x});
SL_ankle_angle(x,:)=SL_ankle_angle(x,:)+[0 0 90];
SL_ankle_angle(x,:)=SL_ankle_angle(x,:).*[1 , -1 , -1];


end

SL_ankle_IntExt=mean(SL_ankle_angle(:,1));
SL_ankle_IE=mean(SL_ankle_angle(:,2));
SL_ankle_FE=mean(SL_ankle_angle(:,3));



%Foot Progression angle

vector=SRC_right_ankle-SRC_right_metatarsal;

for x=1:length(SRC_right_ankle)
    
angulo(x)=Angulos_proy_ortogonal(vector(x,:),SRFP_RefS{x}(1,:), SRFP_RefS{x}(3,:), SRFP_RefS{x}(2,:), [0,0,0])-90;

end

SRFp_ankle_IntExt=mean(angulo);

vector1=SRC_left_ankle-SRC_left_metatarsal;

for x=1:length(SRC_right_ankle)
    
angulo1(x)=Angulos_proy_ortogonal(vector1(x,:),SLFP_RefS{x}(1,:), SLFP_RefS{x}(3,:), SLFP_RefS{x}(2,:), [0,0,0])-90;
angulo1(x)=angulo1(x)*-1;

end

SLFp_ankle_IntExt=mean(angulo1);







%% Ángulo de flexo extension de la rodilla, varo y valgo y rotacion interna y externa

%Utilizando el método de BTS

for x=1:length(SRef_sys_RMuslo_BTS)
    
SR_knee_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_RMuslo_BTS{x}, SRef_sys_RPierna{x});

end

SR_knee_IntExt=mean(SR_knee_angle(:,1));
SR_knee_VV=mean(SR_knee_angle(:,2));
SR_knee_FE=mean(SR_knee_angle(:,3));



for x=1:length(SRef_sys_LMuslo_BTS)
    
SL_knee_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_LMuslo_BTS{x}, SRef_sys_LPierna{x});
SL_knee_angle(x,:)=SL_knee_angle(x,:).*[-1 ,-1, 1];

end

SL_knee_IntExt=mean(SL_knee_angle(:,1));
SL_knee_VV=mean(SL_knee_angle(:,2));
SL_knee_FE=mean(SL_knee_angle(:,3));



%Utilizando el Método de Bell

for x=1:length(SRef_sys_RMuslo_Bell)
    
SR1_knee_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_RMuslo_Bell{x}, SRef_sys_RPierna{x});

end

SR1_knee_IntExt=mean(SR1_knee_angle(:,1));
SR1_knee_VV=mean(SR1_knee_angle(:,2));
SR1_knee_FE=mean(SR1_knee_angle(:,3));



for x=1:length(SRef_sys_LMuslo_Bell)
    
SL1_knee_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_LMuslo_Bell{x}, SRef_sys_LPierna{x});
SL1_knee_angle(x,:)=SL1_knee_angle(x,:).*[-1 ,-1, 1];

end

SL1_knee_IntExt=mean(SL1_knee_angle(:,1));
SL1_knee_VV=mean(SL1_knee_angle(:,2));
SL1_knee_FE=mean(SL1_knee_angle(:,3));



% Método del 25%

for x=1:length(SRef_sys_RMuslo_25)
    
SR2_knee_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_RMuslo_25{x}, SRef_sys_RPierna{x});


end

SR2_knee_IntExt=mean(SR2_knee_angle(:,1));
SR2_knee_VV=mean(SR2_knee_angle(:,2));
SR2_knee_FE=mean(SR2_knee_angle(:,3));



for x=1:length(SRef_sys_LMuslo_25)
    
SL2_knee_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_LMuslo_25{x}, SRef_sys_LPierna{x});
SL2_knee_angle(x,:)=SL2_knee_angle(x,:).*[-1 ,-1, 1];

end

SL2_knee_IntExt=mean(SL2_knee_angle(:,1));
SL2_knee_VV=mean(SL2_knee_angle(:,2));
SL2_knee_FE=mean(SL2_knee_angle(:,3));



%Por método de Davis 

for x=1:length(SRef_sys_RMuslo_Davis)
    
SR3_knee_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_RMuslo_Davis{x}, SRef_sys_RPierna{x});

end

SR3_knee_IntExt=mean(SR3_knee_angle(:,1));
SR3_knee_VV=mean(SR3_knee_angle(:,2));
SR3_knee_FE=mean(SR3_knee_angle(:,3));


for x=1:length(SRef_sys_LMuslo_Davis)
    
SL3_knee_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_LMuslo_Davis{x}, SRef_sys_LPierna{x});
SL3_knee_angle(x,:)=SL3_knee_angle(x,:).*[-1 ,-1, 1];

end

SL3_knee_IntExt=mean(SL3_knee_angle(:,1));
SL3_knee_VV=mean(SL3_knee_angle(:,2));
SL3_knee_FE=mean(SL3_knee_angle(:,3));






%% Calcular el angulo de la cadera

%Con método de BTS
SR_Hip_angle = [];
for x=1:length(SRef_sys_Pelvis)
    
SR_Hip_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, SRef_sys_RMuslo_BTS{x});
SR_Hip_angle(x,:)=SR_Hip_angle(x,:).*[1 ,1 , -1];

end

SR_Hip_IntExt=mean(SR_Hip_angle(:,1));
SR_Hip_AA=mean(SR_Hip_angle(:,2));
SR_Hip_FE=mean(SR_Hip_angle(:,3));



for x=1:length(SRef_sys_Pelvis)
    
SL_Hip_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, SRef_sys_LMuslo_BTS{x});
SL_Hip_angle(x,:)=SL_Hip_angle(x,:).*[-1 ,-1 , -1];


end

SL_Hip_IntExt=mean(SL_Hip_angle(:,1));
SL_Hip_AA=mean(SL_Hip_angle(:,2));
SL_Hip_FE=mean(SL_Hip_angle(:,3));

%Con el método de Bell
 
for x=1:length(SRef_sys_Pelvis)
    
SR1_Hip_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, SRef_sys_RMuslo_Bell{x});
SR1_Hip_angle(x,:)=SR1_Hip_angle(x,:).*[1, 1, -1];

end

SR1_Hip_IntExt=mean(SR1_Hip_angle(:,1));
SR1_Hip_AA=mean(SR1_Hip_angle(:,2));
SR1_Hip_FE=mean(SR1_Hip_angle(:,3));



for x=1:length(SRef_sys_Pelvis)
    
SL1_Hip_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, SRef_sys_LMuslo_Bell{x});
SL1_Hip_angle(x,:)=SL1_Hip_angle(x,:).*[-1, -1, -1];


end

SL1_Hip_IntExt=mean(SL1_Hip_angle(:,1));
SL1_Hip_AA=mean(SL1_Hip_angle(:,2));
SL1_Hip_FE=mean(SL1_Hip_angle(:,3));



%Con el método de 25%

for x=1:length(SRef_sys_Pelvis)
    
SR2_Hip_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, SRef_sys_RMuslo_25{x});
SR2_Hip_angle(x,:)=SR2_Hip_angle(x,:).*[1, 1, -1];

end

SR2_Hip_IntExt=mean(SR2_Hip_angle(:,1));
SR2_Hip_AA=mean(SR2_Hip_angle(:,2));
SR2_Hip_FE=mean(SR2_Hip_angle(:,3));



for x=1:length(SRef_sys_Pelvis)
    
SL2_Hip_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, SRef_sys_LMuslo_25{x});
SL2_Hip_angle(x,:)=SL2_Hip_angle(x,:).*[-1, -1, -1];

end

SL2_Hip_IntExt=mean(SL2_Hip_angle(:,1));
SL2_Hip_AA=mean(SL2_Hip_angle(:,2));
SL2_Hip_FE=mean(SL2_Hip_angle(:,3));




% %Con el método de Davis
 
for x=1:length(SRef_sys_Pelvis)
    
SR3_Hip_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, SRef_sys_RMuslo_Davis{x});
SR3_Hip_angle(x,:)=SR3_Hip_angle(x,:).*[1 ,1 , -1];

end

SR3_Hip_IntExt=mean(SR3_Hip_angle(:,1));
SR3_Hip_AA=mean(SR3_Hip_angle(:,2));
SR3_Hip_FE=mean(SR3_Hip_angle(:,3));



for x=1:length(SRef_sys_Pelvis)
    
SL3_Hip_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, SRef_sys_LMuslo_Davis{x});
SL3_Hip_angle(x,:)=SL3_Hip_angle(x,:).*[-1 ,-1 , -1];


end

SL3_Hip_IntExt=mean(SL3_Hip_angle(:,1));
SL3_Hip_AA=mean(SL3_Hip_angle(:,2));
SL3_Hip_FE=mean(SL3_Hip_angle(:,3));



%% Pelvic Tilt/obliquity and Rotation

for x=1:length(SRef_sys_Pelvis)
    
SL_Pelvis_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Pelvis{x}, RefS_global_marcha);
SL_Pelvis_angle(x,:)=SL_Pelvis_angle(x,:).*[1, -1, -1];

end

SL_Pelvis_Rotation=mean(SL_Pelvis_angle(:,1));
SL_Pelvis_Obliquity=mean(SL_Pelvis_angle(:,2));
SL_Pelvis_Tilt=mean(SL_Pelvis_angle(:,3));



for x=1:length(SRef_sys_Pelvis)
    
SR_Pelvis_angle(x,:)=SL_Pelvis_angle(x,:).*[-1, -1, 1];

end

SR_Pelvis_Rotation=mean(SR_Pelvis_angle(:,1));
SR_Pelvis_Obliquity=mean(SR_Pelvis_angle(:,2));
SR_Pelvis_Tilt=mean(SR_Pelvis_angle(:,3));






%% Trunk Tilt, obliquity and rotation

for x=1:length(SRef_sys_Trunk)
    
SL_Trunk_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Trunk{x}, RefS_global_marcha);
SL_Trunk_angle(x,:)=SL_Trunk_angle(x,:).*[1, -1, -1];


end

SL_Trunk_Rotation=mean(SL_Trunk_angle(:,1));
SL_Trunk_Obliquity=mean(SL_Trunk_angle(:,2));
SL_Trunk_Tilt=mean(SL_Trunk_angle(:,3));



for x=1:length(SRef_sys_Trunk)
    
SR_Trunk_angle(x,:)=SL_Trunk_angle(x,:).*[-1, -1, 1];

end

SR_Trunk_Rotation=mean(SR_Trunk_angle(:,1));
SR_Trunk_Obliquity=mean(SR_Trunk_angle(:,2));
SR_Trunk_Tilt=mean(SR_Trunk_angle(:,3));




%% Spinal flexo-extension, medio-lateral bending and internal external rotation

for x=1:length(SRef_sys_Trunk)
    
SR_Spinal_angle(x,:)=EulerA_ZYX([3 2 1], SRef_sys_Trunk{x}, SRef_sys_Pelvis{x});

end

SR_Spinal_Rotation=mean(SR_Spinal_angle(:,1));
SR_Spinal_Bending=mean(SR_Spinal_angle(:,2));
SR_Spinal_FE=mean(SR_Spinal_angle(:,3));


for x=1:length(SRef_sys_Trunk)
    
SL_Spinal_angle(x,:)=SR_Spinal_angle(x,:).*[-1, -1, 1];

end

SL_Spinal_Rotation=mean(SL_Spinal_angle(:,1));
SL_Spinal_Obliquity=mean(SL_Spinal_angle(:,2));
SL_Spinal_Tilt=mean(SL_Spinal_angle(:,3));

% figure();
% subplot(3,1,1);
% plot(SL_Spinal_angle(:,1));
% subplot(3,1,2);
% plot(SL_Spinal_angle(:,2));
% subplot(3,1,3);
% plot(SL_Spinal_angle(:,3));















%% Walking

%Corrección de los asis por medio de las medidas antropometricas
[R_Asis2,L_Asis2]=WAsis_correction(Names,Dynamic_int_fil);

%Sistemas de referencia global para la marcha

WRefS_global_marcha=WEje_Global_marcha(Names, Dynamic_int_fil, r1_on, r2_on, l1_on, l2_on);

%Events normalized on each gait cycle

enRHSLCycle=Fun_enRHSLCycle(r1_on,r2_on,l1_on,l2_on);
enLHSRCycle=Fun_enLHSRCycle(l1_on,l2_on, r1_on,r2_on);

enRTO=Fun_TO(r1_on,r2_on,r_off);
enLTO=Fun_TO(l1_on,l2_on,l_off);

enLTORCycle=Fun_TO(r1_on,r2_on, l_off);
enRTOLCycle=Fun_TO(l1_on,l2_on, r_off);



%Para el caso dado el ultimo dio mayor al 100%











%% Sistemas de referencia técnicos y centros de rotación de los segmentos

%Rodillas
RC_right_knee=Right_Knee_RC(Names, Dynamic_int_fil, R_epy_distance);
RC_left_knee=Left_Knee_RC(Names, Dynamic_int_fil, L_epy_distance);

%Tobillos
RC_right_ankle=Right_Ankle_RC(Names, Dynamic_int_fil, R_mall_distance);
RC_left_ankle=Left_Ankle_RC(Names, Dynamic_int_fil, L_mall_distance);

%Metatarsos
RC_right_metatarsal=Right_Metatarsal_RC(Names, Dynamic_int_fil, R_mall_distance);
RC_left_metatarsal=Left_Metatarsal_RC(Names, Dynamic_int_fil, L_mall_distance);

%Sistema de referencia de el tronco 

[Ref_sys_Trunk,Trunk_RC]=RefS_Trunk(Names, Dynamic_int_fil);

%Pelvis por los diferentes métodos que existen

%método de BTS

[BTS_RC_right_hip,Eje_Pelvis]=BTS_RPelvis_RC(Names, R_Asis2, L_Asis2, Dynamic_int_fil, RHPTP);
BTS_RC_left_hip=BTS_LPelvis_RC(Names, R_Asis2, L_Asis2, Dynamic_int_fil, LHPTP);

%Por método de Bell

Bell_RC_right_hip=Bell_RPelvis_RC(Names,R_Asis2, L_Asis2,Dynamic_int_fil);
Bell_RC_left_hip=Bell_LPelvis_RC(Names,R_Asis2, L_Asis2, Dynamic_int_fil);

%Con Paper Davis 1991

Leg_length=0.5*(L_real_leg_length+R_real_leg_length);
Davis_RC_right_hip=DavisH_RPelvis_RC(Names, Dynamic_int_fil, R_Asis2, L_Asis2, Leg_length);
Davis_RC_left_hip=DavisH_LPelvis_RC(Names, Dynamic_int_fil, R_Asis2, L_Asis2, Leg_length);

%Método del 25%
%Usando las indicaciones de Songning, el vector entre los dos trocanters y
%moverme sobre este un 25% de la distancia

PC25_RC_right_hip=PerC25_RPelvis_RC(Names, Dynamic_int_fil);
PC25_RC_left_hip=PerC25_LPelvis_RC(Names, Dynamic_int_fil);


%% Reference system foot progression-Walking

for x=1:length(RC_right_ankle)
   
    [RFP_RefS{x},LFP_RefS{x}] =WFP_RefS(WRefS_global_marcha, RC_right_ankle(x,:), RC_right_metatarsal(x,:), RC_left_ankle(x,:), RC_left_metatarsal(x,:));
    
end




%% Centros de masa de segmentos

% Encontrar los centros de masa de los segmentos

%Foot segment

RFoot_CoM=CoM(RC_right_ankle,RC_right_metatarsal,0.44);
LFoot_CoM=CoM(RC_left_ankle,RC_left_metatarsal,0.44);

%Shank segment

RShank_CoM=CoM(RC_right_knee,RC_right_ankle,0.42);
LShank_CoM=CoM(RC_left_knee,RC_left_ankle,0.42);

%Thigh Segment

%With Bell method

RThigh_CoM=CoM(Bell_RC_right_hip, RC_right_knee,0.39);
LThigh_CoM=CoM(Bell_RC_left_hip, RC_left_knee,0.39);

%With 25% method

R2Thigh_CoM=CoM(PC25_RC_right_hip, RC_right_knee,0.39);
L2Thigh_CoM=CoM(PC25_RC_left_hip, RC_left_knee,0.39);

%with helen hayes (davis heel)

R3Thigh_CoM=CoM(Davis_RC_right_hip, RC_right_knee,0.39);
L3Thigh_CoM=CoM(Davis_RC_left_hip, RC_left_knee,0.39);

%With BTS method

R4Thigh_CoM=CoM(BTS_RC_right_hip, RC_right_knee,0.39);
L4Thigh_CoM=CoM(BTS_RC_left_hip, RC_left_knee,0.39);







%% Ejes de referencia anatómicos

%Eje anatómico del Muslo por método de Bell

Ref_sys_RMuslo_Bell=RefS_RMuslo(Names, Dynamic_int_fil, Bell_RC_right_hip, RC_right_knee);

for x=1:length(Ref_sys_RMuslo_Bell)
    X_new = Ref_sys_RMuslo_Bell{x}(1,:);
    Y_new = Ref_sys_RMuslo_Bell{x}(3,:);
    Z_new = Ref_sys_RMuslo_Bell{x}(2,:).*-1;
    Ref_sys_RMuslo_Bell{x}=[X_new;Y_new;Z_new];    
end

Ref_sys_LMuslo_Bell=RefS_LMuslo(Names, Dynamic_int_fil, Bell_RC_left_hip, RC_left_knee);

for x=1:length(Ref_sys_LMuslo_Bell)
    X_new = Ref_sys_LMuslo_Bell{x}(1,:);
    Y_new = Ref_sys_LMuslo_Bell{x}(3,:).*-1;
    Z_new = Ref_sys_LMuslo_Bell{x}(2,:);
    Ref_sys_LMuslo_Bell{x}=[X_new;Y_new;Z_new];    
end

%Eje anatómico del Muslo por método del 25%

Ref_sys_RMuslo_25=RefS_RMuslo(Names, Dynamic_int_fil, PC25_RC_right_hip, RC_right_knee);

for x=1:length(Ref_sys_RMuslo_25)
    X_new = Ref_sys_RMuslo_25{x}(1,:);
    Y_new = Ref_sys_RMuslo_25{x}(3,:);
    Z_new = Ref_sys_RMuslo_25{x}(2,:).*-1;
    Ref_sys_RMuslo_25{x}=[X_new;Y_new;Z_new];    
end

Ref_sys_LMuslo_25=RefS_LMuslo(Names, Dynamic_int_fil, PC25_RC_left_hip, RC_left_knee);

for x=1:length(Ref_sys_LMuslo_25)
    X_new = Ref_sys_LMuslo_25{x}(1,:);
    Y_new = Ref_sys_LMuslo_25{x}(3,:).*-1;
    Z_new = Ref_sys_LMuslo_25{x}(2,:);
    Ref_sys_LMuslo_25{x}=[X_new;Y_new;Z_new];    
end

%Eje anatómico del Muslo por método de Davis 1991

Ref_sys_RMuslo_Davis=RefS_RMuslo(Names, Dynamic_int_fil, Davis_RC_right_hip, RC_right_knee);

for x=1:length(Ref_sys_RMuslo_Davis)
    X_new = Ref_sys_RMuslo_Davis{x}(1,:);
    Y_new = Ref_sys_RMuslo_Davis{x}(3,:);
    Z_new = Ref_sys_RMuslo_Davis{x}(2,:).*-1;
    Ref_sys_RMuslo_Davis{x}=[X_new;Y_new;Z_new];    
end

Ref_sys_LMuslo_Davis=RefS_LMuslo(Names, Dynamic_int_fil, Davis_RC_left_hip, RC_left_knee);

for x=1:length(Ref_sys_LMuslo_Davis)
    X_new = Ref_sys_LMuslo_Davis{x}(1,:);
    Y_new = Ref_sys_LMuslo_Davis{x}(3,:).*-1;
    Z_new = Ref_sys_LMuslo_Davis{x}(2,:);
    Ref_sys_LMuslo_Davis{x}=[X_new;Y_new;Z_new];    
end

%Eje anatómico del Muslo por método de BTS

Ref_sys_RMuslo_BTS=RefS_RMuslo(Names, Dynamic_int_fil, BTS_RC_right_hip, RC_right_knee);

for x=1:length(Ref_sys_RMuslo_BTS)
    X_new = Ref_sys_RMuslo_BTS{x}(1,:);
    Y_new = Ref_sys_RMuslo_BTS{x}(3,:);
    Z_new = Ref_sys_RMuslo_BTS{x}(2,:).*-1;
    Ref_sys_RMuslo_BTS{x}=[X_new;Y_new;Z_new];    
end

Ref_sys_LMuslo_BTS=RefS_LMuslo(Names, Dynamic_int_fil, BTS_RC_left_hip, RC_left_knee);

for x=1:length(Ref_sys_LMuslo_BTS)
    X_new = Ref_sys_LMuslo_BTS{x}(1,:);
    Y_new = Ref_sys_LMuslo_BTS{x}(3,:).*-1;
    Z_new = Ref_sys_LMuslo_BTS{x}(2,:);
    Ref_sys_LMuslo_BTS{x}=[X_new;Y_new;Z_new];    

end

%Eje anatómico de la Pierna

Ref_sys_RPierna=RefS_RPierna(Names, Dynamic_int_fil, RC_right_knee, RC_right_ankle);
Ref_sys_LPierna=RefS_LPierna(Names, Dynamic_int_fil, RC_left_knee, RC_left_ankle);

% Eje anatómico del pie

Ref_sys_RPie_BTS=RefS_RPie_BTS(Names, Dynamic_int_fil, RC_right_ankle, RC_right_metatarsal);
Ref_sys_LPie_BTS=RefS_LPie_BTS(Names, Dynamic_int_fil, RC_left_ankle, RC_left_metatarsal);

%Eje de referencia de la pelvis

Ref_sys_Pelvis = {};
for x=1:length(Eje_Pelvis)
    X_new = Eje_Pelvis{x}(3,:);
    Y_new = Eje_Pelvis{x}(1,:);
    Z_new = Eje_Pelvis{x}(2,:);
    Ref_sys_Pelvis{x}=[X_new;Y_new;Z_new];    
end


%Eje de referencia de el tronco

for x=1:length(Ref_sys_Trunk)
    X_new = Ref_sys_Trunk{x}(3,:);
    Y_new = Ref_sys_Trunk{x}(2,:).*-1;
    Z_new = Ref_sys_Trunk{x}(1,:);
    Ref_sys_Trunk{x}=[X_new;Y_new;Z_new];    
end












 

%% Sacar los ángulos

%Calcular ángulos del pie (dorsi y plantar flexión, inversión-eversión,
%rotación interna y externa)

for x=1:length(Ref_sys_RPierna)
    
R_ankle_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_RPierna{x}, Ref_sys_RPie_BTS{x});
R_ankle_angle(x,:)=R_ankle_angle(x,:)+[0 0 90];
R_ankle_angle(x,:)=R_ankle_angle(x,:).*[-1, 1, -1];

end

R_ankle_IntExt=R_ankle_angle(r1_on:r2_on,1);
R_ankle_IE=R_ankle_angle(r1_on:r2_on,2);
R_ankle_FE=R_ankle_angle(r1_on:r2_on,3);




for x=1:length(Ref_sys_LPierna)
    
L_ankle_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_LPierna{x}, Ref_sys_LPie_BTS{x});
L_ankle_angle(x,:)=L_ankle_angle(x,:)+[0 0 90];
L_ankle_angle(x,:)=L_ankle_angle(x,:).*[1 , -1 , -1];

end
 
L_ankle_IntExt=L_ankle_angle(l1_on:l2_on,1);
L_ankle_IE=L_ankle_angle(l1_on:l2_on,2);
L_ankle_FE=L_ankle_angle(l1_on:l2_on,3);



%%Angulo de progresión del pie 

vector=RC_right_ankle-RC_right_metatarsal;

for x=1:length(RC_right_ankle)
    
angulo(x)=Angulos_proy_ortogonal(vector(x,:),RFP_RefS{x}(1,:), RFP_RefS{x}(3,:), RFP_RefS{x}(2,:), [0,0,0])-90;

end


RFp_ankle_IntExt=angulo(r1_on:r2_on);



vector1=RC_left_ankle-RC_left_metatarsal;

for x=1:length(RC_right_ankle)
    
angulo1(x)=Angulos_proy_ortogonal(vector1(x,:),LFP_RefS{x}(1,:), LFP_RefS{x}(3,:), LFP_RefS{x}(2,:), [0,0,0])-90;
angulo1(x)=angulo1(x)*-1;

end

LFp_ankle_IntExt=angulo1(l1_on:l2_on);



%% Ángulo de flexo extension de la rodilla, varo y valgo y rotacion interna y externa

%Utilizando el método de BTS

for x=1:length(Ref_sys_RMuslo_BTS)
    
R_knee_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_RMuslo_BTS{x}, Ref_sys_RPierna{x});


end

R_knee_IntExt=R_knee_angle(r1_on:r2_on,1);
R_knee_VV=R_knee_angle(r1_on:r2_on,2);
R_knee_FE=R_knee_angle(r1_on:r2_on,3);


for x=1:length(Ref_sys_LMuslo_BTS)
    
L_knee_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_LMuslo_BTS{x}, Ref_sys_LPierna{x});
L_knee_angle(x,:)=L_knee_angle(x,:).*[-1 ,-1, 1];

end

L_knee_IntExt=L_knee_angle(l1_on:l2_on,1);
L_knee_VV=L_knee_angle(l1_on:l2_on,2);
L_knee_FE=L_knee_angle(l1_on:l2_on,3);


%Utilizando el Método de Bell

for x=1:length(Ref_sys_RMuslo_Bell)
    
R1_knee_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_RMuslo_Bell{x}, Ref_sys_RPierna{x});

end

R1_knee_IntExt=R1_knee_angle(r1_on:r2_on,1);
R1_knee_VV=R1_knee_angle(r1_on:r2_on,2);
R1_knee_FE=R1_knee_angle(r1_on:r2_on,3);


for x=1:length(Ref_sys_LMuslo_Bell)
    
L1_knee_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_LMuslo_Bell{x}, Ref_sys_LPierna{x});
L1_knee_angle(x,:)=L1_knee_angle(x,:).*[-1 ,-1, 1];

end

L1_knee_IntExt=L1_knee_angle(l1_on:l2_on,1);
L1_knee_VV=L1_knee_angle(l1_on:l2_on,2);
L1_knee_FE=L1_knee_angle(l1_on:l2_on,3);



% Método del 25%

for x=1:length(Ref_sys_RMuslo_25)
    
R2_knee_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_RMuslo_25{x}, Ref_sys_RPierna{x});

end

R2_knee_IntExt=R2_knee_angle(r1_on:r2_on,1);
R2_knee_VV=R2_knee_angle(r1_on:r2_on,2);
R2_knee_FE=R2_knee_angle(r1_on:r2_on,3);



for x=1:length(Ref_sys_LMuslo_25)
    
L2_knee_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_LMuslo_25{x}, Ref_sys_LPierna{x});
L2_knee_angle(x,:)=L2_knee_angle(x,:).*[-1 ,-1, 1];

end

L2_knee_IntExt=L2_knee_angle(l1_on:l2_on,1);
L2_knee_VV=L2_knee_angle(l1_on:l2_on,2);
L2_knee_FE=L2_knee_angle(l1_on:l2_on,3);



%Por método de Davis 

for x=1:length(Ref_sys_RMuslo_Davis)
    
R3_knee_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_RMuslo_Davis{x}, Ref_sys_RPierna{x});

end

R3_knee_IntExt=R3_knee_angle(r1_on:r2_on,1);
R3_knee_VV=R3_knee_angle(r1_on:r2_on,2);
R3_knee_FE=R3_knee_angle(r1_on:r2_on,3);

for x=1:length(Ref_sys_LMuslo_Davis)
    
L3_knee_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_LMuslo_Davis{x}, Ref_sys_LPierna{x});
L3_knee_angle(x,:)=L3_knee_angle(x,:).*[-1 ,-1, 1];

end

L3_knee_IntExt=L3_knee_angle(l1_on:l2_on,1);
L3_knee_VV=L3_knee_angle(l1_on:l2_on,2);
L3_knee_FE=L3_knee_angle(l1_on:l2_on,3);



%% Calcular el angulo de la cadera

%Con método de BTS

for x=1:length(Ref_sys_Pelvis)
    
R_Hip_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, Ref_sys_RMuslo_BTS{x});
R_Hip_angle(x,:)=R_Hip_angle(x,:).*[1 ,1 , -1];

end

R_hip_IntExt=R_Hip_angle(r1_on:r2_on,1);
R_hip_AA=R_Hip_angle(r1_on:r2_on,2);
R_hip_FE=R_Hip_angle(r1_on:r2_on,3);



for x=1:length(Ref_sys_Pelvis)
    
L_Hip_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, Ref_sys_LMuslo_BTS{x});
L_Hip_angle(x,:)=L_Hip_angle(x,:).*[-1 ,-1 , -1];

end


L_hip_IntExt=L_Hip_angle(l1_on:l2_on,1);
L_hip_AA=L_Hip_angle(l1_on:l2_on,2);
L_hip_FE=L_Hip_angle(l1_on:l2_on,3);


%Con el método de Bell
 
for x=1:length(Ref_sys_Pelvis)
    
R1_Hip_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, Ref_sys_RMuslo_Bell{x});
R1_Hip_angle(x,:)=R1_Hip_angle(x,:).*[1, 1, -1];

end

R1_hip_IntExt=R1_Hip_angle(r1_on:r2_on,1);
R1_hip_AA=R1_Hip_angle(r1_on:r2_on,2);
R1_hip_FE=R1_Hip_angle(r1_on:r2_on,3);



for x=1:length(Ref_sys_Pelvis)
    
L1_Hip_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, Ref_sys_LMuslo_Bell{x});
L1_Hip_angle(x,:)=L1_Hip_angle(x,:).*[-1, -1, -1];

end

L1_hip_IntExt=L1_Hip_angle(l1_on:l2_on,1);
L1_hip_AA=L1_Hip_angle(l1_on:l2_on,2);
L1_hip_FE=L1_Hip_angle(l1_on:l2_on,3);




%Con el método de 25%

for x=1:length(Ref_sys_Pelvis)
    
R2_Hip_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, Ref_sys_RMuslo_25{x});
R2_Hip_angle(x,:)=R2_Hip_angle(x,:).*[1, 1, -1];

end

R2_hip_IntExt=R2_Hip_angle(r1_on:r2_on,1);
R2_hip_AA=R2_Hip_angle(r1_on:r2_on,2);
R2_hip_FE=R2_Hip_angle(r1_on:r2_on,3);


for x=1:length(Ref_sys_Pelvis)
    
L2_Hip_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, Ref_sys_LMuslo_25{x});
L2_Hip_angle(x,:)=L2_Hip_angle(x,:).*[-1, -1, -1];

end

L2_hip_IntExt=L2_Hip_angle(l1_on:l2_on,1);
L2_hip_AA=L2_Hip_angle(l1_on:l2_on,2);
L2_hip_FE=L2_Hip_angle(l1_on:l2_on,3);



% %Con el método de Davis
 
for x=1:length(Ref_sys_Pelvis)
    
R3_Hip_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, Ref_sys_RMuslo_Davis{x});
R3_Hip_angle(x,:)=R3_Hip_angle(x,:).*[1 ,1 , -1];


end

R3_hip_IntExt=R3_Hip_angle(r1_on:r2_on,1);
R3_hip_AA=R3_Hip_angle(r1_on:r2_on,2);
R3_hip_FE=R3_Hip_angle(r1_on:r2_on,3);



for x=1:length(Ref_sys_Pelvis)
    
L3_Hip_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, Ref_sys_LMuslo_Davis{x});
L3_Hip_angle(x,:)=L3_Hip_angle(x,:).*[-1 ,-1 , -1];

end

L3_hip_IntExt=L3_Hip_angle(l1_on:l2_on,1);
L3_hip_AA=L3_Hip_angle(l1_on:l2_on,2);
L3_hip_FE=L3_Hip_angle(l1_on:l2_on,3);




%% Pelvic Tilt/obliquity and Rotation

for x=1:length(Ref_sys_Pelvis)
    
L_Pelvis_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Pelvis{x}, WRefS_global_marcha);
L_Pelvis_angle(x,:)=L_Pelvis_angle(x,:).*[1, -1, -1];

end

L_Pelvis_Rotation=L_Pelvis_angle(l1_on:l2_on,1);
L_Pelvis_Obliquity=L_Pelvis_angle(l1_on:l2_on,2);
L_Pelvis_Tilt=L_Pelvis_angle(l1_on:l2_on,3);



for x=1:length(Ref_sys_Pelvis)
    
R_Pelvis_angle(x,:)=L_Pelvis_angle(x,:).*[-1, -1, 1];

end

R_Pelvis_Rotation=R_Pelvis_angle(r1_on:r2_on,1);
R_Pelvis_Obliquity=R_Pelvis_angle(r1_on:r2_on,2);
R_Pelvis_Tilt=R_Pelvis_angle(r1_on:r2_on,3);



%% Trunk Tilt, obliquity and rotation

for x=1:length(Ref_sys_Trunk)
    
L_Trunk_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Trunk{x}, WRefS_global_marcha);
L_Trunk_angle(x,:)=L_Trunk_angle(x,:).*[1, -1, -1];

end
 
L_Trunk_Rotation=L_Trunk_angle(l1_on:l2_on,1);
L_Trunk_Obliquity=L_Trunk_angle(l1_on:l2_on,2);
L_Trunk_Tilt=L_Trunk_angle(l1_on:l2_on,3);



for x=1:length(Ref_sys_Trunk)
    
R_Trunk_angle(x,:)=L_Trunk_angle(x,:).*[-1, -1, 1];

end

R_Trunk_Rotation=R_Trunk_angle(r1_on:r2_on,1);
R_Trunk_Obliquity=R_Trunk_angle(r1_on:r2_on,2);
R_Trunk_Tilt=R_Trunk_angle(r1_on:r2_on,3);


%% Ángulos del tronco relativos a los ángulos de standing

L_Trunk_Tilt_off=L_Trunk_Tilt-SL_Trunk_Tilt;
L_Trunk_Obliquity_off=L_Trunk_Obliquity-SL_Trunk_Obliquity;
L_Trunk_Rotation_off=L_Trunk_Rotation-SL_Trunk_Rotation;

R_Trunk_Tilt_off=R_Trunk_Tilt-SR_Trunk_Tilt;
R_Trunk_Obliquity_off=R_Trunk_Obliquity-SR_Trunk_Obliquity;
R_Trunk_Rotation_off=R_Trunk_Rotation-SR_Trunk_Rotation;




%% Spinal flexo-extension, medio-lateral bending and internal external rotation

for x=1:length(Ref_sys_Trunk)
    
R_Spinal_angle(x,:)=EulerA_ZYX([3 2 1], Ref_sys_Trunk{x}, Ref_sys_Pelvis{x});

end

R_Spinal_Rotation=R_Spinal_angle(r1_on:r2_on,1);
R_Spinal_Obliquity=R_Spinal_angle(r1_on:r2_on,2);
R_Spinal_Tilt=R_Spinal_angle(r1_on:r2_on,3);



for x=1:length(Ref_sys_Trunk)
    
L_Spinal_angle(x,:)=R_Spinal_angle(x,:).*[-1, -1, 1];

end

L_Spinal_Rotation=L_Spinal_angle(l1_on:l2_on,1);
L_Spinal_Obliquity=L_Spinal_angle(l1_on:l2_on,2);
L_Spinal_Tilt=L_Spinal_angle(l1_on:l2_on,3);

% figure(1);
% subplot(3,1,1);
% plot(L_Spinal_angle(:,1));
% subplot(3,1,2);
% plot(L_Spinal_angle(:,2));
% subplot(3,1,3);
% plot(L_Spinal_angle(:,3));


%% Tiempo de ciclo

tseqRStride=(r2_on-r1_on)/100;
tseqLStride=(l2_on-l1_on)/100;

%% Longitud de paso 

R_stride_length=RStep_length(Names, Dynamic_int_fil, r1_on, r2_on);
L_stride_length=LStep_length(Names, Dynamic_int_fil, l1_on, l2_on);

%% Ancho de paso

m=Find_name(Names,'r heel.X');
n=Find_name(Names,'r heel.Y');
o=Find_name(Names,'r heel.Z');

R_Heel=[Dynamic_int_fil(:,m), Dynamic_int_fil(:,n), Dynamic_int_fil(:,o)];

M=Find_name(Names,'l heel.X');
N=Find_name(Names,'l heel.Y');
O=Find_name(Names,'l heel.Z');

L_Heel=[Dynamic_int_fil(:,M), Dynamic_int_fil(:,N), Dynamic_int_fil(:,O)];


if l1_on>r1_on

dis3=norm(L_Heel(l1_on,:)-R_Heel(r1_on,:));
dis4=norm(L_Heel(l2_on,:)-R_Heel(r2_on,:));
dis1=norm(R_Heel(r2_on,:)-L_Heel(l1_on,:));
dis2=norm(R_Heel(r2_on,:)-L_Heel(l1_on,:));

vector1=(R_Heel(r2_on,:)-L_Heel(l1_on,:))/norm((R_Heel(r2_on,:)-L_Heel(l1_on,:)));
vector2=(R_Heel(r2_on,:)-L_Heel(l1_on,:))/norm((R_Heel(r2_on,:)-L_Heel(l1_on,:)));
vector3=(L_Heel(l1_on,:)-R_Heel(r1_on,:))/norm(L_Heel(l1_on,:)-R_Heel(r1_on,:));
vector4=(L_Heel(l2_on,:)-R_Heel(r2_on,:))/norm(L_Heel(l2_on,:)-R_Heel(r2_on,:));

nRStep=1;
nLStep=2;
tStep=3;

else 

dis1=norm(R_Heel(r1_on,:)-L_Heel(l1_on,:));
dis2=norm(R_Heel(r2_on,:)-L_Heel(l2_on,:));
dis3=norm(L_Heel(l2_on,:)-R_Heel(r1_on,:));
dis4=norm(L_Heel(l2_on,:)-R_Heel(r1_on,:));

vector1=(R_Heel(r1_on,:)-L_Heel(l1_on,:))/norm(R_Heel(r1_on,:)-L_Heel(l1_on,:));
vector2=(R_Heel(r2_on,:)-L_Heel(l2_on,:))/norm(R_Heel(r2_on,:)-L_Heel(l2_on,:));
vector3=(L_Heel(l2_on,:)-R_Heel(r1_on,:))/norm(L_Heel(l2_on,:)-R_Heel(r1_on,:));
vector4=(L_Heel(l2_on,:)-R_Heel(r1_on,:))/norm(L_Heel(l2_on,:)-R_Heel(r1_on,:));

nRStep=2;
nLStep=1;
tStep=3;

%vector 3y4 es el RHStoLHS, los otros dos son LHStoRHS

end


aGAITRHS(1,1)=acosd(dot(vector1,WRefS_global_marcha(2,:)));
aGAITRHS(1,2)=acosd(dot(vector2,WRefS_global_marcha(2,:)));

sinaGAITRHS=[sind(aGAITRHS(1,1)),sind(aGAITRHS(1,2))];
cosaGAITRHS=[cosd(aGAITRHS(1,1)),cosd(aGAITRHS(1,2))];

RStepwidth(1,1)=dis1*sinaGAITRHS(1);
RStepwidth(1,2)=dis2*sinaGAITRHS(2);

aGAITLHS(1,1)=acosd(dot(vector3,WRefS_global_marcha(2,:)));
aGAITLHS(1,2)=acosd(dot(vector4,WRefS_global_marcha(2,:)));

sinaGAITLHS=[sind(aGAITLHS(1,1)),sind(aGAITLHS(1,2))];
cosaGAITLHS=[cosd(aGAITLHS(1,1)),cosd(aGAITLHS(1,2))];

LStepwidth(1,1)=dis3*sinaGAITLHS(1);
LStepwidth(1,2)=dis4*sinaGAITLHS(2);

%Step width averaged in entire cycle

RStepwidth=mean(RStepwidth);
LStepwidth=mean(LStepwidth);

dtotRStepwidth=RStepwidth*nRStep;
dtotLStepwidth=LStepwidth*nLStep;

dtotStepwidth=dtotRStepwidth+dtotLStepwidth;

dmStepwidth=dtotStepwidth/tStep;

%Step calculated at heel strike

dseqRSTep=[dis1*cosaGAITRHS(1) dis2*cosaGAITRHS(2)];
dseqLSTep=[dis3*cosaGAITLHS(1) dis4*cosaGAITLHS(2)];

dmRStep=mean(dseqRSTep);
dmLStrp=mean(dseqLSTep);


%Stance, Swing  and double support phases duration averaged on the entire
%walking sequence

tseqRStance=(r_off-r1_on)/100;
tseqRSwing=(r2_on-r_off)/100;

tseqLStance=(l_off-l1_on)/100;
tseqLSwing=(l2_on-l_off)/100;


if r1_on>l1_on
    
   tseqRDBLStance=(l_off-r1_on)/100;
   tseqLDBLStance=(r_off-l2_on)/100;

else
    
   tseqRDBLStance=(l_off-r2_on)/100;
   tseqLDBLStance=(r_off-l1_on)/100;
   
end

sRStance=tseqRStance/tseqRStride;
sRSwing=tseqRSwing/tseqRStride;
sRDBLStance=tseqRDBLStance/tseqRStride;

sLStance=tseqLStance/tseqLStride;
sLSwing=tseqLSwing/tseqLStride;
sLDBLStance=tseqLDBLStance/tseqLStride;

sRSINGStance=tseqLSwing/tseqRStride;
sLSINGStance=tseqRSwing/tseqLStride;

%Velocidades de paso para pie derecho e izquierdo

RVelocidad=diff(R_Heel);
RVelocidad=RVelocidad.*100;

LVelocidad=diff(L_Heel);
LVelocidad=LVelocidad.*100;


RVelocidad=[RVelocidad(r1_on:r2_on,1),RVelocidad(r1_on:r2_on,2),RVelocidad(r1_on:r2_on,3)];
meanRVelocidad=mean(RVelocidad);

LVelocidad=[LVelocidad(l1_on:l2_on,1),LVelocidad(l1_on:l2_on,2),LVelocidad(l1_on:l2_on,3)];
meanLVelocidad=mean(LVelocidad);

%Stride Velocities
RVel=norm(meanRVelocidad);
LVel=norm(meanLVelocidad);

%Swing velocity
RVSwing=diff(R_Heel);
RVSwing=RVSwing.*100;

RVSwing=[RVSwing(r_off:r2_on,1),RVSwing(r_off:r2_on,2),RVSwing(r_off:r2_on,3)];
mRVSwing=mean(RVSwing);

LVSwing=diff(L_Heel);
LVSwing=LVSwing.*100;

LVSwing=[LVSwing(l_off:l2_on,1),LVSwing(l_off:l2_on,2),LVSwing(l_off:l2_on,3)];
mLVSwing=mean(LVSwing);

RVelSwing=norm(mRVSwing);
LVelSwing=norm(mLVSwing);

nStrides=2;

VtotalStrides=RVel+LVel;
VmTotalStrides=VtotalStrides/nStrides;

%Cadencia

fRcadence=1/tseqRStride*60;
fLcadence=1/tseqLStride*60;

fCadence=fRcadence+fLcadence;
















%%

Final_Data = [];
Final_names = {};
indica = 1;

% Centros de rotacion cadera TODOS LOS METODOS

Final_Data(:,indica) = BTS_RC_right_hip(:,1);
Final_names{indica} = 'BTS_RC_right_hip.X';
indica = indica+1;
Final_Data(:,indica) = BTS_RC_right_hip(:,2);
Final_names{indica} = 'BTS_RC_right_hip.Y';
indica = indica+1;
Final_Data(:,indica) = BTS_RC_right_hip(:,3);
Final_names{indica} = 'BTS_RC_right_hip.Z';
indica = indica+1;

Final_Data(:,indica) = Bell_RC_right_hip(:,1);
Final_names{indica} = 'Bell_RC_right_hip.X';
indica = indica+1;
Final_Data(:,indica) = Bell_RC_right_hip(:,2);
Final_names{indica} = 'Bell_RC_right_hip.Y';
indica = indica+1;
Final_Data(:,indica) = Bell_RC_right_hip(:,3);
Final_names{indica} = 'Bell_RC_right_hip.Z';
indica = indica+1;

Final_Data(:,indica) = Davis_RC_right_hip(:,1);
Final_names{indica} = 'Davis_RC_right_hip.X';
indica = indica+1;
Final_Data(:,indica) = Davis_RC_right_hip(:,2);
Final_names{indica} = 'Davis_RC_right_hip.Y';
indica = indica+1;
Final_Data(:,indica) = Davis_RC_right_hip(:,3);
Final_names{indica} = 'Davis_RC_right_hip.Z';
indica = indica+1;

Final_Data(:,indica) = PC25_RC_right_hip(:,1);
Final_names{indica} = 'PC25_RC_right_hip.X';
indica = indica+1;
Final_Data(:,indica) = PC25_RC_right_hip(:,2);
Final_names{indica} = 'PC25_RC_right_hip.Y';
indica = indica+1;
Final_Data(:,indica) = PC25_RC_right_hip(:,3);
Final_names{indica} = 'PC25_RC_right_hip.Z';
indica = indica+1;


Final_Data(:,indica) = BTS_RC_left_hip(:,1);
Final_names{indica} = 'BTS_RC_left_hip.X';
indica = indica+1;
Final_Data(:,indica) = BTS_RC_left_hip(:,2);
Final_names{indica} = 'BTS_RC_left_hip.Y';
indica = indica+1;
Final_Data(:,indica) = BTS_RC_left_hip(:,3);
Final_names{indica} = 'BTS_RC_left_hip.Z';
indica = indica+1;

Final_Data(:,indica) = Bell_RC_left_hip(:,1);
Final_names{indica} = 'Bell_RC_left_hip.X';
indica = indica+1;
Final_Data(:,indica) = Bell_RC_left_hip(:,2);
Final_names{indica} = 'Bell_RC_left_hip.Y';
indica = indica+1;
Final_Data(:,indica) = Bell_RC_left_hip(:,3);
Final_names{indica} = 'Bell_RC_left_hip.Z';
indica = indica+1;

Final_Data(:,indica) = Davis_RC_left_hip(:,1);
Final_names{indica} = 'Davis_RC_left_hip.X';
indica = indica+1;
Final_Data(:,indica) = Davis_RC_left_hip(:,2);
Final_names{indica} = 'Davis_RC_left_hip.Y';
indica = indica+1;
Final_Data(:,indica) = Davis_RC_left_hip(:,3);
Final_names{indica} = 'Davis_RC_left_hip.Z';
indica = indica+1;

Final_Data(:,indica) = PC25_RC_left_hip(:,1);
Final_names{indica} = 'PC25_RC_left_hip.X';
indica = indica+1;
Final_Data(:,indica) = PC25_RC_left_hip(:,2);
Final_names{indica} = 'PC25_RC_left_hip.Y';
indica = indica+1;
Final_Data(:,indica) = PC25_RC_left_hip(:,3);
Final_names{indica} = 'PC25_RC_left_hip.Z';
indica = indica+1;

% Centros de rotación de la rodilla

Final_Data(:,indica) = RC_right_knee(:,1);
Final_names{indica} = 'RC_right_knee.X';
indica = indica+1;
Final_Data(:,indica) = RC_right_knee(:,2);
Final_names{indica} = 'RC_right_knee.Y';
indica = indica+1;
Final_Data(:,indica) = RC_right_knee(:,3);
Final_names{indica} = 'RC_right_knee.Z';
indica = indica+1;

Final_Data(:,indica) = RC_left_knee(:,1);
Final_names{indica} = 'RC_left_knee.X';
indica = indica+1;
Final_Data(:,indica) = RC_left_knee(:,2);
Final_names{indica} = 'RC_left_knee.Y';
indica = indica+1;
Final_Data(:,indica) = RC_left_knee(:,3);
Final_names{indica} = 'RC_left_knee.Z';
indica = indica+1;


% Centros de rotación del tobillo

Final_Data(:,indica) = RC_right_ankle(:,1);
Final_names{indica} = 'RC_right_ankle.X';
indica = indica+1;
Final_Data(:,indica) = RC_right_ankle(:,2);
Final_names{indica} = 'RC_right_ankle.Y';
indica = indica+1;
Final_Data(:,indica) = RC_right_ankle(:,3);
Final_names{indica} = 'RC_right_ankle.Z';
indica = indica+1;

Final_Data(:,indica) = RC_left_ankle(:,1);
Final_names{indica} = 'RC_left_ankle.X';
indica = indica+1;
Final_Data(:,indica) = RC_left_ankle(:,2);
Final_names{indica} = 'RC_left_ankle.Y';
indica = indica+1;
Final_Data(:,indica) = RC_left_ankle(:,3);
Final_names{indica} = 'RC_left_ankle.Z';
indica = indica+1;

%Centro de rotación metatarsos

Final_Data(:,indica) = RC_right_metatarsal(:,1);
Final_names{indica} = 'RC_right_metatarsal.X';
indica = indica+1;
Final_Data(:,indica) = RC_right_metatarsal(:,2);
Final_names{indica} = 'RC_right_metatarsal.Y';
indica = indica+1;
Final_Data(:,indica) = RC_right_metatarsal(:,3);
Final_names{indica} = 'RC_right_metatarsal.Z';
indica = indica+1;

Final_Data(:,indica) = RC_left_metatarsal(:,1);
Final_names{indica} = 'RC_left_metatarsal.X';
indica = indica+1;
Final_Data(:,indica) = RC_left_metatarsal(:,2);
Final_names{indica} = 'RC_left_metatarsal.Y';
indica = indica+1;
Final_Data(:,indica) = RC_left_metatarsal(:,3);
Final_names{indica} = 'RC_left_metatarsal.Z';
indica = indica+1;

%Centros de masa por todos los metodos del muslo

Final_Data(:,indica) = RThigh_CoM(:,1);
Final_names{indica} = 'RThigh_CoM.X';
indica = indica+1;
Final_Data(:,indica) = RThigh_CoM(:,2);
Final_names{indica} = 'RThigh_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = RThigh_CoM(:,3);
Final_names{indica} = 'RThigh_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = R2Thigh_CoM(:,1);
Final_names{indica} = 'R2Thigh_CoM.X';
indica = indica+1;
Final_Data(:,indica) = R2Thigh_CoM(:,2);
Final_names{indica} = 'R2Thigh_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = R2Thigh_CoM(:,3);
Final_names{indica} = 'R2Thigh_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = R3Thigh_CoM(:,1);
Final_names{indica} = 'R3Thigh_CoM.X';
indica = indica+1;
Final_Data(:,indica) = R3Thigh_CoM(:,2);
Final_names{indica} = 'R3Thigh_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = R3Thigh_CoM(:,3);
Final_names{indica} = 'R3Thigh_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = R4Thigh_CoM(:,1);
Final_names{indica} = 'R4Thigh_CoM.X';
indica = indica+1;
Final_Data(:,indica) = R4Thigh_CoM(:,2);
Final_names{indica} = 'R4Thigh_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = R4Thigh_CoM(:,3);
Final_names{indica} = 'R4Thigh_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = RShank_CoM(:,1);
Final_names{indica} = 'RShank_CoM.X';
indica = indica+1;
Final_Data(:,indica) = RShank_CoM(:,2);
Final_names{indica} = 'RShank_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = RShank_CoM(:,3);
Final_names{indica} = 'RShank_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = RFoot_CoM(:,1);
Final_names{indica} = 'RFoot_CoM.X';
indica = indica+1;
Final_Data(:,indica) = RFoot_CoM(:,2);
Final_names{indica} = 'RFoot_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = RFoot_CoM(:,3);
Final_names{indica} = 'RFoot_CoM.Z';
indica = indica+1;



Final_Data(:,indica) = LThigh_CoM(:,1);
Final_names{indica} = 'LThigh_CoM.X';
indica = indica+1;
Final_Data(:,indica) = LThigh_CoM(:,2);
Final_names{indica} = 'LThigh_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = LThigh_CoM(:,3);
Final_names{indica} = 'LThigh_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = L2Thigh_CoM(:,1);
Final_names{indica} = 'L2Thigh_CoM.X';
indica = indica+1;
Final_Data(:,indica) = L2Thigh_CoM(:,2);
Final_names{indica} = 'L2Thigh_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = L2Thigh_CoM(:,3);
Final_names{indica} = 'L2Thigh_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = L3Thigh_CoM(:,1);
Final_names{indica} = 'L3Thigh_CoM.X';
indica = indica+1;
Final_Data(:,indica) = L3Thigh_CoM(:,2);
Final_names{indica} = 'L3Thigh_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = L3Thigh_CoM(:,3);
Final_names{indica} = 'L3Thigh_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = L4Thigh_CoM(:,1);
Final_names{indica} = 'L4Thigh_CoM.X';
indica = indica+1;
Final_Data(:,indica) = L4Thigh_CoM(:,2);
Final_names{indica} = 'L4Thigh_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = L4Thigh_CoM(:,3);
Final_names{indica} = 'L4Thigh_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = LShank_CoM(:,1);
Final_names{indica} = 'LShank_CoM.X';
indica = indica+1;
Final_Data(:,indica) = LShank_CoM(:,2);
Final_names{indica} = 'LShank_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = LShank_CoM(:,3);
Final_names{indica} = 'LShank_CoM.Z';
indica = indica+1;

Final_Data(:,indica) = LFoot_CoM(:,1);
Final_names{indica} = 'LFoot_CoM.X';
indica = indica+1;
Final_Data(:,indica) = LFoot_CoM(:,2);
Final_names{indica} = 'LFoot_CoM.Y';
indica = indica+1;
Final_Data(:,indica) = LFoot_CoM(:,3);
Final_names{indica} = 'LFoot_CoM.Z';
indica = indica+1;


%Ejes de referencia del muslo

for x = 1:length(Ref_sys_RMuslo_Bell)

Final_Data(x,indica) = Ref_sys_RMuslo_Bell{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_RMuslo_Bell{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_RMuslo_Bell{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_RMuslo_Bell{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_RMuslo_Bell{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_RMuslo_Bell{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_RMuslo_Bell{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_RMuslo_Bell{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_RMuslo_Bell{x}(3,3);

end

Final_names{indica} = 'Ref_sys_RMuslo_Bell_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Bell_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Bell_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RMuslo_Bell_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Bell_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Bell_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RMuslo_Bell_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Bell_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Bell_K.Z';
indica = indica+1;


for x = 1:length(Ref_sys_RMuslo_Davis)

Final_Data(x,indica) = Ref_sys_RMuslo_Davis{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_RMuslo_Davis{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_RMuslo_Davis{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_RMuslo_Davis{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_RMuslo_Davis{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_RMuslo_Davis{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_RMuslo_Davis{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_RMuslo_Davis{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_RMuslo_Davis{x}(3,3);

end

Final_names{indica} = 'Ref_sys_RMuslo_Davis_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Davis_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Davis_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RMuslo_Davis_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Davis_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Davis_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RMuslo_Davis_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Davis_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_Davis_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_RMuslo_25)

Final_Data(x,indica) = Ref_sys_RMuslo_25{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_RMuslo_25{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_RMuslo_25{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_RMuslo_25{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_RMuslo_25{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_RMuslo_25{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_RMuslo_25{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_RMuslo_25{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_RMuslo_25{x}(3,3);

end

Final_names{indica} = 'Ref_sys_RMuslo_25_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_25_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_25_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RMuslo_25_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_25_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_25_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RMuslo_25_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_25_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_25_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_RMuslo_BTS)

Final_Data(x,indica) = Ref_sys_RMuslo_BTS{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_RMuslo_BTS{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_RMuslo_BTS{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_RMuslo_BTS{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_RMuslo_BTS{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_RMuslo_BTS{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_RMuslo_BTS{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_RMuslo_BTS{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_RMuslo_BTS{x}(3,3);

end

Final_names{indica} = 'Ref_sys_RMuslo_BTS_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_BTS_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_BTS_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RMuslo_BTS_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_BTS_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_BTS_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RMuslo_BTS_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_BTS_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RMuslo_BTS_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_LMuslo_Bell)

Final_Data(x,indica) = Ref_sys_LMuslo_Bell{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_LMuslo_Bell{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_LMuslo_Bell{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_LMuslo_Bell{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_LMuslo_Bell{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_LMuslo_Bell{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_LMuslo_Bell{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_LMuslo_Bell{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_LMuslo_Bell{x}(3,3);

end

Final_names{indica} = 'Ref_sys_LMuslo_Bell_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Bell_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Bell_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LMuslo_Bell_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Bell_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Bell_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LMuslo_Bell_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Bell_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Bell_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_LMuslo_Davis)

Final_Data(x,indica) = Ref_sys_LMuslo_Davis{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_LMuslo_Davis{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_LMuslo_Davis{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_LMuslo_Davis{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_LMuslo_Davis{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_LMuslo_Davis{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_LMuslo_Davis{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_LMuslo_Davis{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_LMuslo_Davis{x}(3,3);

end

Final_names{indica} = 'Ref_sys_LMuslo_Davis_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Davis_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Davis_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LMuslo_Davis_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Davis_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Davis_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LMuslo_Davis_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Davis_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_Davis_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_LMuslo_25)

Final_Data(x,indica) = Ref_sys_LMuslo_25{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_LMuslo_25{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_LMuslo_25{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_LMuslo_25{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_LMuslo_25{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_LMuslo_25{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_LMuslo_25{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_LMuslo_25{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_LMuslo_25{x}(3,3);

end

Final_names{indica} = 'Ref_sys_LMuslo_25_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_25_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_25_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LMuslo_25_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_25_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_25_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LMuslo_25_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_25_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_25_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_LMuslo_BTS)

Final_Data(x,indica) = Ref_sys_LMuslo_BTS{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_LMuslo_BTS{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_LMuslo_BTS{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_LMuslo_BTS{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_LMuslo_BTS{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_LMuslo_BTS{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_LMuslo_BTS{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_LMuslo_BTS{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_LMuslo_BTS{x}(3,3);

end

Final_names{indica} = 'Ref_sys_LMuslo_BTS_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_BTS_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_BTS_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LMuslo_BTS_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_BTS_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_BTS_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LMuslo_BTS_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_BTS_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LMuslo_BTS_K.Z';
indica = indica+1;

%Eje de referencia anatómico de la pierna

for x = 1:length(Ref_sys_RPierna)

Final_Data(x,indica) = Ref_sys_RPierna{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_RPierna{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_RPierna{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_RPierna{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_RPierna{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_RPierna{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_RPierna{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_RPierna{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_RPierna{x}(3,3);

end

Final_names{indica} = 'Ref_sys_RPierna_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPierna_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPierna_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RPierna_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPierna_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPierna_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RPierna_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPierna_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPierna_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_LPierna)

Final_Data(x,indica) = Ref_sys_LPierna{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_LPierna{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_LPierna{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_LPierna{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_LPierna{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_LPierna{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_LPierna{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_LPierna{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_LPierna{x}(3,3);

end

Final_names{indica} = 'Ref_sys_LPierna_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPierna_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPierna_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LPierna_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPierna_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPierna_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LPierna_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPierna_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPierna_K.Z';
indica = indica+1;

%Eje de referencia antómico del pie

for x = 1:length(Ref_sys_RPie_BTS)

Final_Data(x,indica) = Ref_sys_RPie_BTS{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_RPie_BTS{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_RPie_BTS{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_RPie_BTS{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_RPie_BTS{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_RPie_BTS{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_RPie_BTS{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_RPie_BTS{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_RPie_BTS{x}(3,3);

end

Final_names{indica} = 'Ref_sys_RPie_BTS_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPie_BTS_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPie_BTS_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RPie_BTS_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPie_BTS_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPie_BTS_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_RPie_BTS_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPie_BTS_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_RPie_BTS_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_LPie_BTS)

Final_Data(x,indica) = Ref_sys_LPie_BTS{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_LPie_BTS{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_LPie_BTS{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_LPie_BTS{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_LPie_BTS{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_LPie_BTS{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_LPie_BTS{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_LPie_BTS{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_LPie_BTS{x}(3,3);

end

Final_names{indica} = 'Ref_sys_LPie_BTS_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPie_BTS_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPie_BTS_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LPie_BTS_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPie_BTS_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPie_BTS_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_LPie_BTS_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPie_BTS_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_LPie_BTS_K.Z';
indica = indica+1;

for x = 1:length(Ref_sys_Pelvis)

Final_Data(x,indica) = Ref_sys_Pelvis{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_Pelvis{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_Pelvis{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_Pelvis{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_Pelvis{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_Pelvis{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_Pelvis{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_Pelvis{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_Pelvis{x}(3,3);

end

Final_names{indica} = 'Ref_sys_Pelvis_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Pelvis_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Pelvis_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_Pelvis_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Pelvis_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Pelvis_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_Pelvis_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Pelvis_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Pelvis_K.Z';
indica = indica+1;


for x = 1:length(Ref_sys_Trunk)

Final_Data(x,indica) = Ref_sys_Trunk{x}(1,1);
Final_Data(x,indica+1) = Ref_sys_Trunk{x}(1,2);
Final_Data(x,indica+2) = Ref_sys_Trunk{x}(1,3);


Final_Data(x,indica+3) = Ref_sys_Trunk{x}(2,1);
Final_Data(x,indica+4) = Ref_sys_Trunk{x}(2,2);
Final_Data(x,indica+5) = Ref_sys_Trunk{x}(2,3);


Final_Data(x,indica+6) = Ref_sys_Trunk{x}(3,1);
Final_Data(x,indica+7) = Ref_sys_Trunk{x}(3,2);
Final_Data(x,indica+8) = Ref_sys_Trunk{x}(3,3);

end

Final_names{indica} = 'Ref_sys_Trunk_I.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Trunk_I.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Trunk_I.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_Trunk_J.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Trunk_J.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Trunk_J.Z';
indica = indica+1;

Final_names{indica} = 'Ref_sys_Trunk_K.X';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Trunk_K.Y';
indica = indica+1;
Final_names{indica} = 'Ref_sys_Trunk_K.Z';
indica = indica+1;

% % Angulos del pie
% 
% Final_Data(:,indica) = R_ankle_angle(:,1);
% Final_names{indica} = 'R_ankle_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = R_ankle_angle(:,2);
% Final_names{indica} = 'R_ankle_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = R_ankle_angle(:,3);
% Final_names{indica} = 'R_ankle_angle.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L_ankle_angle(:,1);
% Final_names{indica} = 'L_ankle_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = L_ankle_angle(:,2);
% Final_names{indica} = 'L_ankle_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = L_ankle_angle(:,3);
% Final_names{indica} = 'L_ankle_angle.Z';
% indica = indica+1;
% 
% %Angulo progesion del pie
% 
% Final_Data(:,indica) = angulo(1,:)';
% Final_names{indica} = 'R_progression_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = angulo(1,:)';
% Final_names{indica} = 'R_progression_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = angulo(1,:)';
% Final_names{indica} = 'R_progression_angle.Z';
% indica = indica+1;
% 
% Final_Data(:,indica) = angulo1(1,:)';
% Final_names{indica} = 'L_progression_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = angulo1(1,:)';
% Final_names{indica} = 'L_progression_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = angulo1(1,:)';
% Final_names{indica} = 'L_progression_angle.Z';
% indica = indica+1;
% 
% 
% %Ángulos de la rodilla método de BTS
% 
% Final_Data(:,indica) = R_knee_angle(:,1);
% Final_names{indica} = 'R_knee_angle_BTS.X';
% indica = indica+1;
% Final_Data(:,indica) = R_knee_angle(:,2);
% Final_names{indica} = 'R_knee_angle_BTS.Y';
% indica = indica+1;
% Final_Data(:,indica) = R_knee_angle(:,3);
% Final_names{indica} = 'R_knee_angle_BTS.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L_knee_angle(:,1);
% Final_names{indica} = 'L_knee_angle_BTS.X';
% indica = indica+1;
% Final_Data(:,indica) = L_knee_angle(:,2);
% Final_names{indica} = 'L_knee_angle_BTS.Y';
% indica = indica+1;
% Final_Data(:,indica) = L_knee_angle(:,3);
% Final_names{indica} = 'L_knee_angle_BTS.Z';
% indica = indica+1;
% 
% %ängulos de la Rodilla por metodo de Bell
% 
% Final_Data(:,indica) = R1_knee_angle(:,1);
% Final_names{indica} = 'R_knee_angle_Bell.X';
% indica = indica+1;
% Final_Data(:,indica) = R1_knee_angle(:,2);
% Final_names{indica} = 'R_knee_angle_Bell.Y';
% indica = indica+1;
% Final_Data(:,indica) = R1_knee_angle(:,3);
% Final_names{indica} = 'R_knee_angle_Bell.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L1_knee_angle(:,1);
% Final_names{indica} = 'L_knee_angle_Bell.X';
% indica = indica+1;
% Final_Data(:,indica) = L1_knee_angle(:,2);
% Final_names{indica} = 'L_knee_angle_Bell.Y';
% indica = indica+1;
% Final_Data(:,indica) = L1_knee_angle(:,3);
% Final_names{indica} = 'L_knee_angle_Bell.Z';
% indica = indica+1;
% 
%  %ängulos de la Rodilla por metodo de Davis
% 
% Final_Data(:,indica) = R3_knee_angle(:,1);
% Final_names{indica} = 'R_knee_angle_Davis.X';
% indica = indica+1;
% Final_Data(:,indica) = R3_knee_angle(:,2);
% Final_names{indica} = 'R_knee_angle_Davis.Y';
% indica = indica+1;
% Final_Data(:,indica) = R3_knee_angle(:,3);
% Final_names{indica} = 'R_knee_angle_Davis.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L3_knee_angle(:,1);
% Final_names{indica} = 'L_knee_angle_Davis.X';
% indica = indica+1;
% Final_Data(:,indica) = L3_knee_angle(:,2);
% Final_names{indica} = 'L_knee_angle_Davis.Y';
% indica = indica+1;
% Final_Data(:,indica) = L3_knee_angle(:,3);
% Final_names{indica} = 'L_knee_angle_Davis.Z';
% indica = indica+1;   
% 
%  %ängulos de la Rodilla por metodo de 25%
% 
% Final_Data(:,indica) = R2_knee_angle(:,1);
% Final_names{indica} = 'R_knee_angle_25.X';
% indica = indica+1;
% Final_Data(:,indica) = R2_knee_angle(:,2);
% Final_names{indica} = 'R_knee_angle_25.Y';
% indica = indica+1;
% Final_Data(:,indica) = R2_knee_angle(:,3);
% Final_names{indica} = 'R_knee_angle_25.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L2_knee_angle(:,1);
% Final_names{indica} = 'L_knee_angle_25.X';
% indica = indica+1;
% Final_Data(:,indica) = L2_knee_angle(:,2);
% Final_names{indica} = 'L_knee_angle_25.Y';
% indica = indica+1;
% Final_Data(:,indica) = L2_knee_angle(:,3);
% Final_names{indica} = 'L_knee_angle_25.Z';
% indica = indica+1; 
% 
% %Angulo cadera por metodo de BTS
% 
% Final_Data(:,indica) = R_Hip_angle(:,1);
% Final_names{indica} = 'R_Hip_angle_BTS.X';
% indica = indica+1;
% Final_Data(:,indica) = R_Hip_angle(:,2);
% Final_names{indica} = 'R_Hip_angle_BTS.Y';
% indica = indica+1;
% Final_Data(:,indica) = R_Hip_angle(:,3);
% Final_names{indica} = 'R_Hip_angle_BTS.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L_Hip_angle(:,1);
% Final_names{indica} = 'L_Hip_angle_BTS.X';
% indica = indica+1;
% Final_Data(:,indica) = L_Hip_angle(:,2);
% Final_names{indica} = 'L_Hip_angle_BTS.Y';
% indica = indica+1;
% Final_Data(:,indica) = L_Hip_angle(:,3);
% Final_names{indica} = 'L_Hip_angle_BTS.Z';
% indica = indica+1; 
% 
% %Angulo cadera por metodo de Bell
% 
% Final_Data(:,indica) = R1_Hip_angle(:,1);
% Final_names{indica} = 'R_Hip_angle_Bell.X';
% indica = indica+1;
% Final_Data(:,indica) = R1_Hip_angle(:,2);
% Final_names{indica} = 'R_Hip_angle_Bell.Y';
% indica = indica+1;
% Final_Data(:,indica) = R1_Hip_angle(:,3);
% Final_names{indica} = 'R_Hip_angle_Bell.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L1_Hip_angle(:,1);
% Final_names{indica} = 'L_Hip_angle_Bell.X';
% indica = indica+1;
% Final_Data(:,indica) = L1_Hip_angle(:,2);
% Final_names{indica} = 'L_Hip_angle_Bell.Y';
% indica = indica+1;
% Final_Data(:,indica) = L1_Hip_angle(:,3);
% Final_names{indica} = 'L_Hip_angle_Bell.Z';
% indica = indica+1;
% 
% %Angulo cadera por metodo de Davis
% 
% Final_Data(:,indica) = R3_Hip_angle(:,1);
% Final_names{indica} = 'R_Hip_angle_Davis.X';
% indica = indica+1;
% Final_Data(:,indica) = R3_Hip_angle(:,2);
% Final_names{indica} = 'R_Hip_angle_Davis.Y';
% indica = indica+1;
% Final_Data(:,indica) = R3_Hip_angle(:,3);
% Final_names{indica} = 'R_Hip_angle_Davis.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L3_Hip_angle(:,1);
% Final_names{indica} = 'L_Hip_angle_Davis.X';
% indica = indica+1;
% Final_Data(:,indica) = L3_Hip_angle(:,2);
% Final_names{indica} = 'L_Hip_angle_Davis.Y';
% indica = indica+1;
% Final_Data(:,indica) = L3_Hip_angle(:,3);
% Final_names{indica} = 'L_Hip_angle_Davis.Z';
% indica = indica+1;
% 
% %Angulo cadera por metodo de 25%
% 
% Final_Data(:,indica) = R2_Hip_angle(:,1);
% Final_names{indica} = 'R_Hip_angle_25.X';
% indica = indica+1;
% Final_Data(:,indica) = R2_Hip_angle(:,2);
% Final_names{indica} = 'R_Hip_angle_25.Y';
% indica = indica+1;
% Final_Data(:,indica) = R2_Hip_angle(:,3);
% Final_names{indica} = 'R_Hip_angle_25.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L2_Hip_angle(:,1);
% Final_names{indica} = 'L_Hip_angle_25.X';
% indica = indica+1;
% Final_Data(:,indica) = L2_Hip_angle(:,2);
% Final_names{indica} = 'L_Hip_angle_25.Y';
% indica = indica+1;
% Final_Data(:,indica) = L2_Hip_angle(:,3);
% Final_names{indica} = 'L_Hip_angle_25.Z';
% indica = indica+1;
% 
% %Angulo de la pelvis
% 
% Final_Data(:,indica) = R_Pelvis_angle(:,1);
% Final_names{indica} = 'R_Pelvis_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = R_Pelvis_angle(:,2);
% Final_names{indica} = 'R_Pelvis_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = R_Pelvis_angle(:,3);
% Final_names{indica} = 'R_Pelvis_angle.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L_Pelvis_angle(:,1);
% Final_names{indica} = 'L_Pelvis_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = L_Pelvis_angle(:,2);
% Final_names{indica} = 'L_Pelvis_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = L_Pelvis_angle(:,3);
% Final_names{indica} = 'L_Pelvis_angle.Z';
% indica = indica+1;
% 
% %Angulos del tronco
% 
% Final_Data(:,indica) = R_Trunk_angle(:,1);
% Final_names{indica} = 'R_Trunk_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = R_Trunk_angle(:,2);
% Final_names{indica} = 'R_Trunk_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = R_Trunk_angle(:,3);
% Final_names{indica} = 'R_Trunk_angle.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L_Trunk_angle(:,1);
% Final_names{indica} = 'L_Trunk_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = L_Trunk_angle(:,2);
% Final_names{indica} = 'L_Trunk_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = L_Trunk_angle(:,3);
% Final_names{indica} = 'L_Trunk_angle.Z';
% indica = indica+1;
% 
% %Angulo de la espina
% 
% %Angulos del tronco
% 
% Final_Data(:,indica) = R_Spinal_angle(:,1);
% Final_names{indica} = 'R_Spinal_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = R_Spinal_angle(:,2);
% Final_names{indica} = 'R_Spinal_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = R_Spinal_angle(:,3);
% Final_names{indica} = 'R_Spinal_angle.Z';
% indica = indica+1;
% 
% 
% Final_Data(:,indica) = L_Spinal_angle(:,1);
% Final_names{indica} = 'L_Spinal_angle.X';
% indica = indica+1;
% Final_Data(:,indica) = L_Spinal_angle(:,2);
% Final_names{indica} = 'L_Spinal_angle.Y';
% indica = indica+1;
% Final_Data(:,indica) = L_Spinal_angle(:,3);
% Final_names{indica} = 'L_Spinal_angle.Z';
% indica = indica+1;



FData{Pac}=Final_Data;
FNames{Pac}=Final_names;


end

% figure();
% subplot(3,1,1);
% plot(L_Spinal_angle(:,1));
% subplot(3,1,2);
% plot(L_Spinal_angle(:,2));
% subplot(3,1,3);
% plot(L_Spinal_angle(:,3));



% for x=1:length(jRF)
%    
%     X(x,:)=jRF{x}(1,:);    
%        
% end
% 
% figure();
% subplot(3,1,1);
% plot(X(:,1));
% subplot(3,1,2);
% plot(X(:,2));
% subplot(3,1,3);
% plot(X(:,3));

 %  figure();
%  plot(L_knee_IntExt)



