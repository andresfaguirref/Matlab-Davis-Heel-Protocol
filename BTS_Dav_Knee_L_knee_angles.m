function L_knee_Angles = BTS_Dav_Knee_L_knee_angles(data,Ref_L_thigh, Ref_L_calf,p)

    for x=1:length(Ref_L_thigh)
        
        Angles = Euler_A([3 2 1], Ref_L_thigh{x}, Ref_L_calf{x});
        
        L_knee_Angles(x,1) = Angles(1) * -1;
        L_knee_Angles(x,2) = Angles(2) * -1;
        L_knee_Angles(x,3) = Angles(3);
        
    end
    
    if p == 1
        
       figure()
       subplot(3,1,1)
       plot(data(:,2),L_knee_Angles(:,1))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Knee Fle-Ext (X)')
       subplot(3,1,2)
       plot(data(:,2),L_knee_Angles(:,2))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Knee Abd-Add (Y)')
       subplot(3,1,3)
       plot(data(:,2),L_knee_Angles(:,3))
       ylabel('Measure (m)')
       xlabel('time (s)')
       title('Left Knee Rot Ext-Int (Z)')
       
    end

end