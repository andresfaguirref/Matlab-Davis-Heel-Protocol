function t = Point_to_plane(punto, vector_normal, punto_plano)

d=-1*(vector_normal(1)*punto_plano(1)+vector_normal(2)*punto_plano(2)+vector_normal(3)*punto_plano(3));
t=-1*((d+dot(vector_normal,punto))/(vector_normal(1)^2+vector_normal(2)^2+vector_normal(3)^2));

end