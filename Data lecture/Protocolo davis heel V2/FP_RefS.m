function [eje1,eje2] =FP_RefS( rgait, rc_ankle1, rc_metatarso1, rc_ankle2, rc_metatarso2)

vector1=rc_ankle1-rc_metatarso1;
vector1=vector1/norm(vector1);

nuevo_punto1=Proy_punto_plano(vector1, rgait(3,:), [0 0 0]);

%nuevo_punto1=nuevo_punto1/norm(nuevo_punto1);

vRAU1=cross(nuevo_punto1, rgait(3,:));
vRAU1=vRAU1/norm(vRAU1);

vRAV1=cross(rgait(3,:),vRAU1);
vRAV1=vRAV1/norm(vRAV1);

rgait(3,:)=rgait(3,:)/norm(rgait(3,:));


vector2=rc_ankle2-rc_metatarso2;
vector2=vector2/norm(vector2);

nuevo_punto2=Proy_punto_plano(vector2, rgait(3,:), [0 0 0]);
%nuevo_punto2=nuevo_punto2/norm(nuevo_punto2);

vRAU2=cross(nuevo_punto2, rgait(3,:));
vRAU2=vRAU2/norm(vRAU2);

vRAV2=cross(rgait(3,:),vRAU2);
vRAV2=vRAV2/norm(vRAV2);

rgait(3,:)=rgait(3,:)/norm(rgait(3,:));


eje1=[-vRAU1; -vRAV1; rgait(3,:)];
eje2=[-vRAU2; -vRAV2; rgait(3,:)];

end