function L_thigh_CG = BTS_Dav_Heel_L_thigh_CG(L_hip, L_knee,time,p)

    Vec = L_knee - L_hip;
    
    estimation_Thigh_CG = 0.455;
    
    for x=1:length(L_hip(:,1))  
 
        distance = norm(Vec(x,:));
        L_thigh_CG(x,:) = move_point_with_dis_vec(L_hip(x,:), distance*estimation_Thigh_CG, Vec(x,:)); 
   
    end
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(time,L_thigh_CG(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left thigh gravity center X')
       subplot(3,1,2)
       plot(time,L_thigh_CG(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left thigh gravity center Y')
       subplot(3,1,3)
       plot(time,L_thigh_CG(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left thigh gravity center Z')
       
    end

end