function R_thigh_RefSys = BTS_Dav_Heel_R_thigh_RefSys(data, names ,R_Thigh_CG,R_hip_RC,R_knee_RC,p);

    R_thigh_RefSys = {};
    
    posX = Find_name(names,'r bar 1.X');
    posY = Find_name(names,'r bar 1.Y');
    posZ = Find_name(names,'r bar 1.Z');
    
    R_bar1 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    plot_i = [];
    plot_j = [];
    plot_k = [];
    for x=1:length(R_Thigh_CG)
        
        Ref_sys2 = Ref_sys_with_3points(R_knee_RC(x,:),R_hip_RC(x,:),R_bar1(x,:));
        
        Ref_sys = Rot_RefSys_one_axis(Ref_sys2, 90, 'X', 0);
        
        R_thigh_RefSys{x}(1,:) = Ref_sys(1,:);
        plot_i(x,:) = Ref_sys(1,:);
        
        
        R_thigh_RefSys{x}(2,:) = Ref_sys(2,:);
        plot_j(x,:) = Ref_sys(2,:);
        
        R_thigh_RefSys{x}(3,:) = Ref_sys(3,:);
        plot_k(x,:) = Ref_sys(3,:);
        
        R_thigh_RefSys{x}(4,:) = R_Thigh_CG(x,:);
        
    end

    if p == 1
        
       figure()
       subplot(3,3,1)
       plot(data(:,2),plot_i(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys I X')
       subplot(3,3,4)
       plot(data(:,2),plot_i(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys I Y')
       subplot(3,3,7)
       plot(data(:,2),plot_i(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys I Z')
       
       subplot(3,3,2)
       plot(data(:,2),plot_j(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys J X')
       subplot(3,3,5)
       plot(data(:,2),plot_j(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys J Y')
       subplot(3,3,8)
       plot(data(:,2),plot_j(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys J Z')
       
       subplot(3,3,3)
       plot(data(:,2),plot_k(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys K X')
       subplot(3,3,6)
       plot(data(:,2),plot_k(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys K Y')
       subplot(3,3,9)
       plot(data(:,2),plot_k(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right thigh RefSys K Z')
       
    end
    
end