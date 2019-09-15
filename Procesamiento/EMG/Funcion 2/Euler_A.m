function Euler_Ang = Euler_A(ord, ref_g, ref_l)
    
    st = ord(1);
    nd = ord(2);
    rd = ord(3);
    
    Linea_nodos = cross(ref_g(st,:),ref_l(rd,:));
    Linea_nodos = Linea_nodos/norm(Linea_nodos);

        
    verificar1 = acosd(dot(Linea_nodos, ref_g(rd,:)));
    
    if verificar1 < 90
        e1 = -1*acosd(dot(Linea_nodos,ref_g(nd,:)));
    else 
        e1 = acosd(dot(Linea_nodos,ref_g(nd,:)));
    end
    
       
    verificar2 = Point_to_plane(ref_l(rd,:), ref_g(st,:), [0 0 0]);
    
    if verificar2 <0
        e2 = -(90 - acosd(dot(ref_g(st,:),ref_l(rd,:))));
    else
        e2 = (90 - acosd(dot(-ref_g(st,:),ref_l(rd,:))));
    end
    
    verificar3 = Point_to_plane(ref_l(nd,:), ref_g(st,:), [0 0 0]);
    
    if verificar3 < 0
         e3 = acosd(dot(Linea_nodos,ref_l(nd,:)));
    else 
         e3 = -acosd(dot(Linea_nodos,ref_l(nd,:)));
    end
    
    Euler_Ang = [e3 e2 e1];
    
end