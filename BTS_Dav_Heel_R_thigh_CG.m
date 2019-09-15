function R_thigh_CG = BTS_Dav_Heel_R_thigh_CG(R_hip, R_knee,time,p)

    Vec = R_knee - R_hip;
    
    estimation_Thigh_CG = 0.455;
    
    for x=1:length(R_hip(:,1))  
 
        distance = norm(Vec(x,:));
        R_thigh_CG(x,:) = move_point_with_dis_vec(R_hip(x,:), distance*estimation_Thigh_CG, Vec(x,:)); 
   
    end
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(time,R_thigh_CG(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh gravity center X')
       subplot(3,1,2)
       plot(time,R_thigh_CG(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh gravity center Y')
       subplot(3,1,3)
       plot(time,R_thigh_CG(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh gravity center Z')
       
    end

end