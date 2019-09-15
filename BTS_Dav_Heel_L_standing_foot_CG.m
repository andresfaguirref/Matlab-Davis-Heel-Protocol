function L_foot_CG = BTS_Dav_Heel_L_standing_foot_CG(data, names ,S_L_metatarsus_RC,p)

    posX = Find_name(names,'l heel.X');
    posY = Find_name(names,'l heel.Y');
    posZ = Find_name(names,'l heel.Z');
    
    L_heel = [data(:,posX) data(:,posY) data(:,posZ)];

    Vec = S_L_metatarsus_RC - L_heel;
    
    estimation_Foot_CG = 0.442;
    
    for x=1:length(L_heel(:,1))  
 
        distance = norm(Vec(x,:));
        L_foot_CG(x,:) = move_point_with_dis_vec(L_heel(x,:), distance*estimation_Foot_CG, Vec(x,:)); 
   
    end
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(data(:,2),L_foot_CG(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot gravity center X')
       subplot(3,1,2)
       plot(data(:,2),L_foot_CG(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot gravity center Y')
       subplot(3,1,3)
       plot(data(:,2),L_foot_CG(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot gravity center Z')
       
    end

end