function eje_ref=Eje_referencia(punto1,punto2,punto3)

Vector_x=punto2-punto1;
Vector_int=punto3-punto1;

Vector_z=cross(Vector_x,Vector_int);
Vector_y=cross(Vector_z, Vector_x);

eje_ref=[Vector_x/norm(Vector_x); Vector_y/norm(Vector_y);Vector_z/norm(Vector_z)];

end