function R_Ankle_Angles = BTS_Dav_Heel_R_ankle_angles(data,Ref_R_calf, Ref_R_foot,p)

    for x=1:length(Ref_R_calf)
        
        Angles = Euler_A([3 2 1], Ref_R_calf{x}, Ref_R_foot{x});
        
        R_Ankle_Angles(x,1) = Angles(1) * -1;
        R_Ankle_Angles(x,2) = Angles(2);
        R_Ankle_Angles(x,3) = (Angles(3) + 90) * -1;
        
    end
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(data(:,2),R_Ankle_Angles(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right Ankle Fle-Ext (X)')
       subplot(3,1,2)
       plot(data(:,2),R_Ankle_Angles(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right Ankle Abd-Add (Y)')
       subplot(3,1,3)
       plot(data(:,2),R_Ankle_Angles(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Right Ankle Rot Ext-Int (Z)')
       
    end

end