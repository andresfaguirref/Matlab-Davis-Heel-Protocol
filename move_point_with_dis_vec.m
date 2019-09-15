function New_point = move_point_with_dis_vec(initial_point,distance,unit_vector)
    
    magnitude = norm(unit_vector);
    
    unit_vector = unit_vector/magnitude;
    
    Xx =  distance * unit_vector(1) + initial_point(1);
    Yy =  distance * unit_vector(2) + initial_point(2);
    Zz =  distance * unit_vector(3) + initial_point(3);
    
    New_point = [Xx Yy Zz];

end
