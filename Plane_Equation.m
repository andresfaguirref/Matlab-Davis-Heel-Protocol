function Plane_Equ = Plane_Equation(vector_nor, plane_ponit)

vector_nor=vector_nor/norm(vector_nor);

d=-1*(vector_nor(1)*plane_ponit(1)+vector_nor(2)*plane_ponit(2)+vector_nor(3)*plane_ponit(3));

Plane_Equ = [vector_nor(1) vector_nor(2) vector_nor(3) d]

end