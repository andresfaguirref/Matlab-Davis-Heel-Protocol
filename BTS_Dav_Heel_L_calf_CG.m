function L_calf_CG = BTS_Dav_Heel_L_calf_CG(L_knee, L_ankle,time,p)

    Vec = L_ankle - L_knee;
    
    estimation_Calf_CG = 0.405;
    
    for x=1:length(L_knee(:,1))  
 
        distance = norm(Vec(x,:));
        L_calf_CG(x,:) = move_point_with_dis_vec(L_knee(x,:), distance*estimation_Calf_CG, Vec(x,:)); 
   
    end
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(time,L_calf_CG(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left calf gravity center X')
       subplot(3,1,2)
       plot(time,L_calf_CG(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left calf gravity center Y')
       subplot(3,1,3)
       plot(time,L_calf_CG(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left calf gravity center Z')
       
    end

end