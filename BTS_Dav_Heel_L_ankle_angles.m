function L_Ankle_Angles = BTS_Dav_Heel_L_ankle_angles(data,Ref_L_calf, Ref_L_foot,p)

    for x=1:length(Ref_L_calf)
        
        Angles = Euler_A([3 2 1], Ref_L_calf{x}, Ref_L_foot{x});
        
         L_Ankle_Angles(x,1) = Angles(1);
         L_Ankle_Angles(x,2) = Angles(2) * -1;
         L_Ankle_Angles(x,3) = (Angles(3) + 90) * -1;
      
    end
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(data(:,2),L_Ankle_Angles(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Ankle Fle-Ext (X)')
       subplot(3,1,2)
       plot(data(:,2),L_Ankle_Angles(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Ankle Abd-Add (Y)')
       subplot(3,1,3)
       plot(data(:,2),L_Ankle_Angles(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Ankle Rot Ext-Int (Z)')
       
    end

end