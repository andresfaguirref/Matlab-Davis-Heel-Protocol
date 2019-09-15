a = load('Walking3Markers.mat');
data = a.Walking3Markers;

a =  load('Walking3names.mat');
data_names = a.Walking3names;

X_postions = Find_all_names(data_names,'.X');
Y_postions = Find_all_names(data_names,'.Y');
Z_postions = Find_all_names(data_names,'.Z');

point = data(:,24:26);
point2 = data(:,48:50);

for x=1:length(X_postions)
    X_data(:,x) = data(:,X_postions(x));
    Y_data(:,x) = data(:,Y_postions(x));
    Z_data(:,x) = data(:,Z_postions(x));
end
%%

max_X = max(max(X_data));
min_X = min(min(X_data));

max_Y = max(max(Y_data));
min_Y = min(min(Y_data));

max_Z = max(max(Z_data));
min_Z = min(min(Z_data));

figure()
ax = axes('XLim',[min_X-0.1*abs(min_X) max_X*1.1],'YLim',[min_Z-0.1*abs(min_Z) max_Z*1.1],'ZLim',[min_Y-0.1*abs(min_Y) max_Y*1.1]);
%ax = axes('XLim',[-1.5 3],'YLim',[-1 1],'ZLim',[0 1]);

view(3)
grid on
axis equal;

hold on
xlabel('X')
ylabel('Z')
zlabel('Y')

%[x,y,z] = sphere(20);

[x,y,z] = cylinder([2 2]);

rad = 0.05;
Initial_position = [0 0 0];

h = surface(x*rad,y*rad,z*rad,'FaceColor','red');
t1 = hgtransform('Parent',ax);
t2 = hgtransform('Parent',ax);
set(h,'Parent',t1)

h2 = copyobj(h,t2);

% [x,y,z] = cylinder([-0.2 2]);
% c = surface(x*rad,y*rad,z*rad,'FaceColor','red');
% t3 = hgtransform('Parent',ax);
% set(c,'Parent',t3)

Txy = makehgtform('translate',[-1.5 -1.5 0]);
set(t2,'Matrix',Txy)

[n,m] = size(point);


h3 = line('XData',[1 2],'YData',[1 2],'ZData',[1 2]);
%h3 = surface(mmm,mmm,mmm,'FaceColor','red');
t3 = hgtransform('Parent',ax);
set(h3,'Parent',t3)

 for i = 1:n
    
    %%xxx = t*i/100;
    %%yyy = t;
    %%zzz = t;
    mmm = line('XData',[i i+1],'YData',[i i+1],'ZData',[i i+1]);
    set(t3,'Parent',mmm)
    %translation = makehgtform('translate',Xx(i),Yy(i),Zz(i));
    translation1 = makehgtform('translate',point(i,1),point(i,3),point(i,2));
    set(t1,'Matrix',translation1)
    %drawnow

   
    translation2 = makehgtform('translate',point2(i,1),point2(i,3),point2(i,2));
    set(t2,'Matrix',translation2)
    %plot3(xxx,yyy,zzz,'Color','red');
    
    drawnow
    
    
    
    %pause(1/100)
end


