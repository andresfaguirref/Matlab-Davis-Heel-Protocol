function L_metatarsus_Rot_Cen = BTS_Dav_Heel_L_metatarsus_RC(data,names,malleolus_distance,p)

    posX = Find_name(names,'l knee 2.X');
    posY = Find_name(names,'l knee 2.Y');
    posZ = Find_name(names,'l knee 2.Z');
    
    L_Knee2 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'l bar 2.X');
    posY = Find_name(names,'l bar 2.Y');
    posZ = Find_name(names,'l bar 2.Z');
    
    L_bar2 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'l mall.X');
    posY = Find_name(names,'l mall.Y');
    posZ = Find_name(names,'l mall.Z');
    
    L_mall = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'l met.X');
    posY = Find_name(names,'l met.Y');
    posZ = Find_name(names,'l met.Z');
    
    L_met = [data(:,posX) data(:,posY) data(:,posZ)];
    
    [n m] = size(L_mall);
    
    L_metatarsus_Rot_Cen1 = [];
    for x=1:n
        
        Ref_Ankle_sys = Ref_sys_with_3points(L_mall(x,:),L_Knee2(x,:),L_bar2(x,:));
        L_metatarsus_Rot_Cen1(x,:) = move_point_with_dis_vec(L_met(x,:),malleolus_distance/-2,Ref_Ankle_sys(2,:));
        
    end
    L_metatarsus_Rot_Cen = L_metatarsus_Rot_Cen1;
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(data(:,2),L_metatarsus_Rot_Cen1(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Metatarsus center rotation X')
       subplot(3,1,2)
       plot(data(:,2),L_metatarsus_Rot_Cen1(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Metatarsus center rotation Y')
       subplot(3,1,3)
       plot(data(:,2),L_metatarsus_Rot_Cen1(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Metatarsus center rotation Z')
    end

end