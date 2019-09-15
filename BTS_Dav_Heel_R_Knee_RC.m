function R_Knee_Rot_Cen = BTS_Dav_Heel_R_Knee_RC(data,names,condyles_distance,p)
    
    posX = Find_name(names,'r knee 1.X');
    posY = Find_name(names,'r knee 1.Y');
    posZ = Find_name(names,'r knee 1.Z');
    
    R_Knee1 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'r bar 1.X');
    posY = Find_name(names,'r bar 1.Y');
    posZ = Find_name(names,'r bar 1.Z');
    
    R_bar1 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'r thigh.X');
    posY = Find_name(names,'r thigh.Y');
    posZ = Find_name(names,'r thigh.Z');
    
    R_thigh = [data(:,posX) data(:,posY) data(:,posZ)];
    
    [n m] = size(R_thigh);
    
    R_Knee_Rot_Cen1 = [];
    for x=1:n
        
        Ref_Knee_sys = Ref_sys_with_3points(R_Knee1(x,:),R_thigh(x,:),R_bar1(x,:));
        R_Knee_Rot_Cen1(x,:) = move_point_with_dis_vec(R_Knee1(x,:),condyles_distance/-2,Ref_Knee_sys(2,:));
        
    end
    R_Knee_Rot_Cen = R_Knee_Rot_Cen1;
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(data(:,2),R_Knee_Rot_Cen1(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right knee center rotation X')
       subplot(3,1,2)
       plot(data(:,2),R_Knee_Rot_Cen1(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right knee center rotation Y')
       subplot(3,1,3)
       plot(data(:,2),R_Knee_Rot_Cen1(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right knee center rotation Z')
    end
    
end