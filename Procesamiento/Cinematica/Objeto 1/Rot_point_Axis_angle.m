function Roted_point = Rot_point_Axis_angle(axis,deg_angle,point)

angle = deg_angle*pi/180;

if angle > 0
    Roted_point = dot(axis,point)*axis + cos(angle)*(point - dot(axis,point)*axis) + sin(angle)*cross(axis,point); 
else
    Roted_point = dot(axis,point)*axis + cos(angle)*(point - dot(axis,point)*axis) - sin(angle)*cross(axis,point); 
end

end