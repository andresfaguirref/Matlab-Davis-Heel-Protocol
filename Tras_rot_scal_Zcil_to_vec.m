function [translation rotation scaling] = Tras_rot_scal_Zcil_to_vec(vec,point_ini)

    vec_n = vec/(norm(vec));
          
    node = cross([0 0 1],vec);
    node_n = node/norm(node);
    angle1 = acos(dot(node_n,[1 0 0]));
          
    angle2 = acos(dot(vec_n,[0 0 1]));
    anglev = acos(dot(node_n,[0 1 0]));
    if (anglev >= pi/2)
        angle1 = (-1)*angle1;
    end
        
    translation = makehgtform('translate',point_ini(1)...
                ,point_ini(2),point_ini(3));
    zrotation = makehgtform('zrotate',angle1);
    xrotation = makehgtform('xrotate',angle2);
    rotation = zrotation*xrotation;
    scaling = makehgtform('scale',norm(vec));
          
end