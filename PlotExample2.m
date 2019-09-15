clf;
figure(1)
rotate3d on
ax = axes('XLim',[0,50],'YLim',[-20,20],'ZLim',[-20,20]);
grid on
axis equal;
view(3);
hold on
xlabel('x')
ylabel('y')
zlabel('z')
%the plane
patch([15 -15 -15 15], [15 15 -15 -15], [0 0 0 0],'red');
hold on
 % the arc
  r=4;
  teta=-pi:0.01:0;
  xn=r*cos(teta);
  yn=r*sin(teta);
  h(1)=plot3(zeros(1,numel(xn)),xn-0.5,yn+5,'green');
  hold on
% the spheres
[xs ,ys, zs]=sphere;
h(2)=surf(xs,ys+3,zs+5,'FaceColor','green');
hold on
h(3)=surf(xs,ys-3,zs+5,'FaceColor','blue');
hold on
t1 = hgtransform('Parent',ax);
t2 = hgtransform('Parent',ax);
set(h,'Parent',t1);
drawnow
% rotate sphere1 holding the arc
    % Form z-axis rotation matrix
    Rz = makehgtform('zrotate',r1);
    % Set transforms for both transform objects
    set(t1,'Matrix',Rz)
    drawnow