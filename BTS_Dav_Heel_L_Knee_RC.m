function L_Knee_Rot_Cen = BTS_Dav_Heel_L_Knee_RC(data,names,condyles_distance,p)
    
    posX = Find_name(names,'l knee 1.X');
    posY = Find_name(names,'l knee 1.Y');
    posZ = Find_name(names,'l knee 1.Z');
    
    L_Knee = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'l bar 1.X');
    posY = Find_name(names,'l bar 1.Y');
    posZ = Find_name(names,'l bar 1.Z');
    
    L_bar1 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'l thigh.X');
    posY = Find_name(names,'l thigh.Y');
    posZ = Find_name(names,'l thigh.Z');
    
    L_thigh = [data(:,posX) data(:,posY) data(:,posZ)];
    
    [n m] = size(L_thigh);
    
    L_Knee_Rot_Cen1 = [];
    for x=1:n
        
        Ref_Knee_sys = Ref_sys_with_3points(L_Knee(x,:),L_thigh(x,:),L_bar1(x,:));
        L_Knee_Rot_Cen1(x,:) = move_point_with_dis_vec(L_Knee(x,:),condyles_distance/-2,Ref_Knee_sys(2,:));
        
    end
    L_Knee_Rot_Cen = L_Knee_Rot_Cen1;
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(data(:,2),L_Knee_Rot_Cen1(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left knee center rotation X')
       subplot(3,1,2)
       plot(data(:,2),L_Knee_Rot_Cen1(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left knee center rotation Y')
       subplot(3,1,3)
       plot(data(:,2),L_Knee_Rot_Cen1(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left knee center rotation Z')
    end
    
end