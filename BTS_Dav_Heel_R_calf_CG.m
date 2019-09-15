function R_calf_CG = BTS_Dav_Heel_R_calf_CG(R_knee, R_ankle,time,p)

    Vec = R_ankle - R_knee;
    
    estimation_Calf_CG = 0.405;
    
    for x=1:length(R_knee(:,1))  
 
        distance = norm(Vec(x,:));
        R_calf_CG(x,:) = move_point_with_dis_vec(R_knee(x,:), distance*estimation_Calf_CG, Vec(x,:)); 
   
    end
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(time,R_calf_CG(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right calf gravity center X')
       subplot(3,1,2)
       plot(time,R_calf_CG(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right calf gravity center Y')
       subplot(3,1,3)
       plot(time,R_calf_CG(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right calf gravity center Z')
       
    end

end