function global_point = Local_to_global_with_point_refsys(ref_point,ref_sys,local_point)

    New_X_point = move_point_with_dis_vec(ref_point,local_point(1),ref_sys(1,:));
    
    New_Y_point = move_point_with_dis_vec(New_X_point,local_point(2),ref_sys(2,:));
    
    global_point = move_point_with_dis_vec(New_Y_point,local_point(3),ref_sys(3,:));
    
end