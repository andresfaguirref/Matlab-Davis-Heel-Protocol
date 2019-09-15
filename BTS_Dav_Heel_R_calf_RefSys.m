function R_calf_RefSys = BTS_Dav_Heel_R_calf_RefSys(data, names ,R_Calf_CG,R_Knee_RC,R_Ankle_RC,p);

    R_calf_RefSys = {};
    
    posX = Find_name(names,'r bar 2.X');
    posY = Find_name(names,'r bar 2.Y');
    posZ = Find_name(names,'r bar 2.Z');
    
    R_bar2 = [data(:,posX) data(:,posY) data(:,posZ)];
    
    plot_i = [];
    plot_j = [];
    plot_k = [];
    for x=1:length(R_Calf_CG)
        
        Ref_sys2 = Ref_sys_with_3points(R_Ankle_RC(x,:),R_Knee_RC(x,:),R_bar2(x,:));
        
        Ref_sys = Rot_RefSys_one_axis(Ref_sys2, 90, 'X', 0);
        
        R_calf_RefSys{x}(1,:) = Ref_sys(1,:);
        plot_i(x,:) = Ref_sys(1,:);
        
        
        R_calf_RefSys{x}(2,:) = Ref_sys(2,:);
        plot_j(x,:) = Ref_sys(2,:);
        R_calf_RefSys{x}(3,:) = Ref_sys(3,:);
        plot_k(x,:) = Ref_sys(3,:);
        
        R_calf_RefSys{x}(4,:) = R_Calf_CG(x,:);
        
    end

    if p == 1
        
       figure()
       subplot(3,3,1)
       plot(data(:,2),plot_i(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys I X')
       subplot(3,3,4)
       plot(data(:,2),plot_i(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys I Y')
       subplot(3,3,7)
       plot(data(:,2),plot_i(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys I Z')
       
       subplot(3,3,2)
       plot(data(:,2),plot_j(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys J X')
       subplot(3,3,5)
       plot(data(:,2),plot_j(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys J Y')
       subplot(3,3,8)
       plot(data(:,2),plot_j(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys J Z')
       
       subplot(3,3,3)
       plot(data(:,2),plot_k(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys K X')
       subplot(3,3,6)
       plot(data(:,2),plot_k(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys K Y')
       subplot(3,3,9)
       plot(data(:,2),plot_k(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right carlf RefSys K Z')
       
    end
    
end