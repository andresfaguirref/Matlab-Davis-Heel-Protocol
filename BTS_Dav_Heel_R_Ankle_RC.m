function R_Ankle_Rot_Cen = BTS_Dav_Heel_R_Ankle_RC(data,names,malleolus_distance,p)

    posX = Find_name(names,'r knee 2.X');
    posY = Find_name(names,'r knee 2.Y');
    posZ = Find_name(names,'r knee 2.Z');
    
    R_Knee2 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'r bar 2.X');
    posY = Find_name(names,'r bar 2.Y');
    posZ = Find_name(names,'r bar 2.Z');
    
    R_bar2 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'r mall.X');
    posY = Find_name(names,'r mall.Y');
    posZ = Find_name(names,'r mall.Z');
    
    R_mall = [data(:,posX) data(:,posY) data(:,posZ)];
    
    [n m] = size(R_mall);
    
    R_Ankle_Rot_Cen1 = [];
    for x=1:n
        
        Ref_Ankle_sys = Ref_sys_with_3points(R_mall(x,:),R_Knee2(x,:),R_bar2(x,:));
        R_Ankle_Rot_Cen1(x,:) = move_point_with_dis_vec(R_mall(x,:),malleolus_distance/-2,Ref_Ankle_sys(2,:));
        
    end
    R_Ankle_Rot_Cen = R_Ankle_Rot_Cen1;
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(data(:,2),R_Ankle_Rot_Cen1(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right Ankle center rotation X')
       subplot(3,1,2)
       plot(data(:,2),R_Ankle_Rot_Cen1(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right Ankle center rotation Y')
       subplot(3,1,3)
       plot(data(:,2),R_Ankle_Rot_Cen1(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right Ankle center rotation Z')
    end

end