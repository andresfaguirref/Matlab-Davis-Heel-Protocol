function RefSys_rot = Rot_RefSys_one_axis(RefSys, angle_degrees, axis, p)

 RefSys_rot = [];
 [n m] = size(RefSys);
 
 if (axis == 'x' || axis == 'X' || axis == 1)
    
     RefSys_rot(1,:) = RefSys(1,:);
     RefSys_rot(2,:) = Rot_point_Axis_angle(RefSys(1,:),angle_degrees,RefSys(2,:));
     RefSys_rot(3,:) = cross(RefSys_rot(1,:),RefSys_rot(2,:));
     
     if n > 3
        RefSys_rot(4,:) = RefSys(4,:);
     end
 
 elseif (axis == 'y' || axis == 'Y' || axis == 2)
     
     
     RefSys_rot(1,:) = Rot_point_Axis_angle(RefSys(2,:),angle_degrees,RefSys(1,:));
     RefSys_rot(2,:) = RefSys(2,:);
     RefSys_rot(3,:) = cross(RefSys_rot(1,:),RefSys_rot(2,:));
     
     if n > 3
        RefSys_rot(4,:) = RefSys(4,:);
     end
     
 elseif (axis == 'z' || axis == 'Z' || axis == 3)
     
     RefSys_rot(1,:) = Rot_point_Axis_angle(RefSys(3,:),angle_degrees,RefSys(1,:));
     RefSys_rot(3,:) = RefSys(3,:);
     RefSys_rot(2,:) = cross(RefSys_rot(3,:),RefSys_rot(1,:));
     
     if n > 3
        RefSys_rot(4,:) = RefSys(4,:);
     end
 
 else
    fprintf('The axis " %s" is not an aviable imput', axis); 
    return
 end

      if p == 1
        
        R1 = RefSys;
        R2 = RefSys_rot;
 
        part = 100;

        t = R1(1,1)/part;
        x1 = [0:t:R1(1,1)];
        t = R1(1,2)/part;
        y1 = [0:t:R1(1,2)];
        t = R1(1,3)/part;
        z1 = [0:t:R1(1,3)];

        t = R1(2,1)/part;
        x2 = [0:t:R1(2,1)];
        t = R1(2,2)/part;
        y2 = [0:t:R1(2,2)];
        t = R1(2,3)/part;
        z2 = [0:t:R1(2,3)];

        t = R1(3,1)/part;
        x3 = [0:t:R1(3,1)];
        t = R1(3,2)/part;
        y3 = [0:t:R1(3,2)];
        t = R1(3,3)/part;
        z3 = [0:t:R1(3,3)];

        g = [0:1/part:1];
        no =  g - g;

        t = R2(1,1)/part;
        x11 = [0:t:R2(1,1)];
        t = R2(1,2)/part;
        y11 = [0:t:R2(1,2)];
        t = R2(1,3)/part;
        z11 = [0:t:R2(1,3)];

        t = R2(2,1)/part;
        x22 = [0:t:R2(2,1)];
        t = R2(2,2)/part;
        y22 = [0:t:R2(2,2)];
        t = R2(2,3)/part;
        z22 = [0:t:R2(2,3)];

        t = R2(3,1)/part;
        x33 = [0:t:R2(3,1)];
        t = R2(3,2)/part;
        y33 = [0:t:R2(3,2)];
        t = R2(3,3)/part;
        z33 = [0:t:R2(3,3)];

        ln = cross(R2(1,:),R1(3,:));

        t = 2*ln(1)/part;
        nx = [-ln(1):t:ln(1)];
        t = 2*ln(2)/part;
        ny = [-ln(2):t:ln(2)];
        t = 2*ln(3)/part;
        nz = [-ln(3):t:ln(3)]; 


        t = 1/part;
        xo = [0:t:1];
        
        if x1 == 0; x1 = no; end
        if y1 == 0; y1 = no; end
        if z1 == 0; z1 = no; end
        
        if x2 == 0; x2 = no; end
        if y2 == 0; y2 = no; end
        if z2 == 0; z2 = no; end
        
        if x3 == 0; x3 = no; end
        if y3 == 0; y3 = no; end
        if z3 == 0; z3 = no; end
        
        if x11 == 0; x11 = no; end
        if y11 == 0; y11 = no; end
        if z11 == 0; z11 = no; end
        
        if x22 == 0; x22 = no; end
        if y22 == 0; y22 = no; end
        if z22 == 0; z22 = no; end
        
        if x33 == 0; x33 = no; end
        if y33 == 0; y33 = no; end
        if z33 == 0; z33 = no; end
        
        figure()
        scatter3(x1,y1,z1,'r>')
        
        text(x1(part), y1(part)*1.2, z1(part)*1.2,'X_1');
        
        ylabel('Y Axis')
        xlabel('X Axis')
        zlabel('Z Axis')
        grid('off')
        hold on
        scatter3(x2,y2,z2,'b>')
        text(x2(part), y2(part)*1.1, z2(part)*1.1,'Y_1');
        scatter3(x3,y3,z3,'g>')
        text(x3(part), y3(part)*1.1, z3(part)*1.1,'Z_1');
        
        scatter3(x11,y11,z11,'r.')
        text(x11(part), y11(part), z11(part)*1.1,'X_2');
        scatter3(x22,y22,z22,'b.')
        text(x22(part), y22(part), z22(part)*1.1,'Y_2')
        scatter3(x33,y33,z33,'g.')
        text(x33(part), y33(part), z33(part)*1.1,'Z_2')

        scatter3(xo,no,no,'k')
        text(1.1, 0, 0,'X_global')
        scatter3(no,xo,no,'k')
        text(0, 1.1, 0,'Y_global')
        scatter3(no,no,xo,'k')
        text(0, 0, 1.1,'Z_global')

        %scatter3(nx,ny,nz,'y')
        %text(nx(part), ny(part), nz(part),'N_line')
         
     end    
 
 
end