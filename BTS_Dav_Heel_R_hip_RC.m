function [R_hip , L_hip] = BTS_Dav_Heel_R_hip_RC(dinamic, names, rasis2, lasis2, Local_R_hip, Local_L_hip,p)

%Me toca hallar el punto medio entre las espinas ilíacas con la distancia y
%mover punto

R_Asis=rasis2;

L_Asis=lasis2;

X=Find_name(names,'sacrum.X');
Y=Find_name(names,'sacrum.Y');
Z=Find_name(names,'sacrum.Z');

Sacrum=[dinamic(:,X), dinamic(:,Y), dinamic(:,Z)];

punto_nuevo=[];
vector=L_Asis-R_Asis;


for x=1:length(L_Asis(:,1))  
 
    distancia=norm(vector(x,:));
   
    punto_nuevo(x,:) = move_point_with_dis_vec(R_Asis(x,:), distancia/2, vector(x,:)); 
    
    
end

punto_mueve_x=[]; 
punto_mueve_y=[];
punto_mueve_z=[];

for x=1:length(L_Asis(:,1))
  
eje_ref = Ref_sys_with_3points(Sacrum(x,:),punto_nuevo(x,:),L_Asis(x,:));
    
punto_mueve_x(x,:) = move_point_with_dis_vec(punto_nuevo(x,:), Local_R_hip(1),eje_ref(1,:)); 
punto_mueve_y(x,:) = move_point_with_dis_vec(punto_mueve_x(x,:), Local_R_hip(2),eje_ref(2,:));
punto_mueve_z(x,:) = move_point_with_dis_vec(punto_mueve_y(x,:), Local_R_hip(3),eje_ref(3,:));
    
end

R_hip = punto_mueve_z;

punto_mueve_x=[]; 
punto_mueve_y=[];
punto_mueve_z=[];

for x=1:length(L_Asis(:,1))
    
punto_mueve_x(x,:) = move_point_with_dis_vec(punto_nuevo(x,:), Local_L_hip(1),eje_ref(1,:)); 
punto_mueve_y(x,:) = move_point_with_dis_vec(punto_mueve_x(x,:), Local_L_hip(2),eje_ref(2,:));
punto_mueve_z(x,:) = move_point_with_dis_vec(punto_mueve_y(x,:), Local_L_hip(3),eje_ref(3,:));
    
end

L_hip = punto_mueve_z;

 if p == 1
        
       figure()
       subplot(3,1,1)
       plot(dinamic(:,2),R_hip(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right hip center rotation X')
       subplot(3,1,2)
       plot(dinamic(:,2),R_hip(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right hip center rotation Y')
       subplot(3,1,3)
       plot(dinamic(:,2),R_hip(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right hip center rotation Z')
       
       
       figure()
       subplot(3,1,1)
       plot(dinamic(:,2),L_hip(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left hip center rotation X')
       subplot(3,1,2)
       plot(dinamic(:,2),L_hip(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left hip center rotation Y')
       subplot(3,1,3)
       plot(dinamic(:,2),L_hip(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left hip center rotation Z')
    end



end