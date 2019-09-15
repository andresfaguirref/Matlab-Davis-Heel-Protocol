function L_Foot_RefSys = BTS_Dav_Heel_L_foot_RefSys(data, names, L_Foot_CG, L_Ankle_RC, L_metatarsus_RC, p)

    L_Foot_RefSys = {};
    
    posX = Find_name(names,'l heel.X');
    posY = Find_name(names,'l heel.Y');
    posZ = Find_name(names,'l heel.Z');
    
    L_heel = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'l mall.X');
    posY = Find_name(names,'l mall.Y');
    posZ = Find_name(names,'l mall.Z');
    
    L_mall = [data(:,posX) data(:,posY) data(:,posZ)];
    
    plot_i = [];
    plot_j = [];
    plot_k = [];
    
    for x=1:length(L_Foot_CG)
        
        Ref_sys(1,:) =  L_heel(x,:) - L_metatarsus_RC(x,:);
        Ref_sys(1,:) = Ref_sys(1,:)/norm(Ref_sys(1,:));
        
        mid_vec = L_mall(x,:) - L_Ankle_RC(x,:);
        
        Ref_sys(2,:) = cross(mid_vec,Ref_sys(1,:));
        Ref_sys(2,:) = Ref_sys(2,:)/norm(Ref_sys(2,:));
        
        Ref_sys(3,:) = cross(Ref_sys(1,:),Ref_sys(2,:));
        Ref_sys(3,:) = Ref_sys(3,:)/norm(Ref_sys(3,:));
        
        L_Foot_RefSys{x}(1,:) = Ref_sys(1,:);
        plot_i(x,:) = Ref_sys(1,:);
        
        L_Foot_RefSys{x}(2,:) = Ref_sys(2,:);
        plot_j(x,:) = Ref_sys(2,:);
        
        L_Foot_RefSys{x}(3,:) = Ref_sys(3,:);
        plot_k(x,:) = Ref_sys(3,:);
        
        L_Foot_RefSys{x}(4,:) = L_Foot_CG(x,:);
        
    end

    if p == 1
        
       figure()
       subplot(3,3,1)
       plot(data(:,2),plot_i(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys I X')
       subplot(3,3,4)
       plot(data(:,2),plot_i(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys I Y')
       subplot(3,3,7)
       plot(data(:,2),plot_i(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys I Z')
       
       subplot(3,3,2)
       plot(data(:,2),plot_j(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys J X')
       subplot(3,3,5)
       plot(data(:,2),plot_j(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys J Y')
       subplot(3,3,8)
       plot(data(:,2),plot_j(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys J Z')
       
       subplot(3,3,3)
       plot(data(:,2),plot_k(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys K X')
       subplot(3,3,6)
       plot(data(:,2),plot_k(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys K Y')
       subplot(3,3,9)
       plot(data(:,2),plot_k(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left foot RefSys K Z')
       
    end
    
end