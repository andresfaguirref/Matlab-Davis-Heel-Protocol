
part = 100;

t = R1(1,1)/part;
x1 = [0:t:R1(1,1)];
t = R1(1,2)/part;
y1 = [0:t:R1(1,2)];
t = R1(1,3)/part;
z1 = [0:t:R1(1,3)];

t = R1(2,1)/part;
x2 = [0:t:R1(2,1)];
t = R1(2,2)/part;
y2 = [0:t:R1(2,2)];
t = R1(2,3)/part;
z2 = [0:t:R1(2,3)];


t = R1(3,1)/part;
x3 = [0:t:R1(3,1)];
t = R1(3,2)/part;
y3 = [0:t:R1(3,2)];
t = R1(3,3)/part;
z3 = [0:t:R1(3,3)];

g = [0:1/part:1];
no =  g - g;

t = R2(1,1)/part;
x11 = [0:t:R2(1,1)];
t = R2(1,2)/part;
y11 = [0:t:R2(1,2)];
t = R2(1,3)/part;
z11 = [0:t:R2(1,3)];

t = R2(2,1)/part;
x22 = [0:t:R2(2,1)];
t = R2(2,2)/part;
y22 = [0:t:R2(2,2)];
t = R2(2,3)/part;
z22 = [0:t:R2(2,3)];

t = R2(3,1)/part;
x33 = [0:t:R2(3,1)];
t = R2(3,2)/part;
y33 = [0:t:R2(3,2)];
t = R2(3,3)/part;
z33 = [0:t:R2(3,3)];

ln = cross(R2(1,:),R1(3,:));

t = 2*ln(1)/part;
nx = [-ln(1):t:ln(1)];
t = 2*ln(2)/part;
ny = [-ln(2):t:ln(2)];
t = 2*ln(3)/part;
nz = [-ln(3):t:ln(3)]; 


t = 1/part;
xo = [0:t:1];



figure()
scatter3(x1,y1,z1,'r>')
hold on
scatter3(x2,y2,z2,'b>')
scatter3(x3,y3,z3,'g>')

scatter3(x11,y11,z11,'r.')
scatter3(x22,y22,z22,'b.')
scatter3(x33,y33,z33,'g.')

scatter3(xo,no,no,'k')
scatter3(no,xo,no,'k')
scatter3(no,no,xo,'k')

scatter3(nx,ny,no,'y')


