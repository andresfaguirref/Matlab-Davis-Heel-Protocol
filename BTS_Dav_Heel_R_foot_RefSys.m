function R_Foot_RefSys = BTS_Dav_Heel_R_foot_RefSys(data, names, R_Foot_CG, R_Ankle_RC, R_metatarsus_RC, p)

    R_Foot_RefSys = {};
    
    posX = Find_name(names,'r heel.X');
    posY = Find_name(names,'r heel.Y');
    posZ = Find_name(names,'r heel.Z');
    
    R_heel = [data(:,posX) data(:,posY) data(:,posZ)];
    
    posX = Find_name(names,'r mall.X');
    posY = Find_name(names,'r mall.Y');
    posZ = Find_name(names,'r mall.Z');
    
    R_mall = [data(:,posX) data(:,posY) data(:,posZ)];
    
    plot_i = [];
    plot_j = [];
    plot_k = [];
    
    for x=1:length(R_Foot_CG)
        
        Ref_sys(1,:) =  R_heel(x,:) - R_metatarsus_RC(x,:);
        Ref_sys(1,:) = Ref_sys(1,:)/norm(Ref_sys(1,:));
        
        mid_vec = R_Ankle_RC(x,:) - R_mall(x,:);
        
        Ref_sys(2,:) = cross(mid_vec,Ref_sys(1,:));
        Ref_sys(2,:) = Ref_sys(2,:)/norm(Ref_sys(2,:));
        
        Ref_sys(3,:) = cross(Ref_sys(1,:),Ref_sys(2,:));
        Ref_sys(3,:) = Ref_sys(3,:)/norm(Ref_sys(3,:));
        
        R_Foot_RefSys{x}(1,:) = Ref_sys(1,:);
        plot_i(x,:) = Ref_sys(1,:);
        
        R_Foot_RefSys{x}(2,:) = Ref_sys(2,:);
        plot_j(x,:) = Ref_sys(2,:);
        
        R_Foot_RefSys{x}(3,:) = Ref_sys(3,:);
        plot_k(x,:) = Ref_sys(3,:);
        
        R_Foot_RefSys{x}(4,:) = R_Foot_CG(x,:);
        
    end

    if p == 1
        
       figure()
       subplot(3,3,1)
       plot(data(:,2),plot_i(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys I X')
       subplot(3,3,4)
       plot(data(:,2),plot_i(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys I Y')
       subplot(3,3,7)
       plot(data(:,2),plot_i(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys I Z')
       
       subplot(3,3,2)
       plot(data(:,2),plot_j(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys J X')
       subplot(3,3,5)
       plot(data(:,2),plot_j(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys J Y')
       subplot(3,3,8)
       plot(data(:,2),plot_j(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys J Z')
       
       subplot(3,3,3)
       plot(data(:,2),plot_k(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys K X')
       subplot(3,3,6)
       plot(data(:,2),plot_k(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys K Y')
       subplot(3,3,9)
       plot(data(:,2),plot_k(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Rigth foot RefSys K Z')
       
    end
    
end