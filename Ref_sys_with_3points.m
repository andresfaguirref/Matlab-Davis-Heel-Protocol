function Ref_sys = Ref_sys_with_3points(p1,p2,p3)

    X_axis = p2 - p1;
    X_axis = X_axis / norm(X_axis);
    
    vec_int = p3 - p1;
    Z_axis = cross(X_axis,vec_int);
    
    Z_axis = Z_axis / norm(Z_axis);
    
    Y_axis = cross(Z_axis,X_axis);
    
    Y_axis = Y_axis / norm(Y_axis);
    
    a = [X_axis Y_axis Z_axis];
    
    Ref_sys = vec2mat(a,3);  
    
end