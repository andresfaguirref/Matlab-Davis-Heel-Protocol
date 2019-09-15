t = 0:pi/50:10*pi;
    x = sin(t);
    y = cos(t);
    z = t;
    ah = axes;
    set(ah,'XLim',[min(x) max(x)],'YLim',[min(y) max(y)],...
        'ZLim',[min(z) max(z)]);
    plot3(x,y,z,'Color','red');
    hold on;
    view(3);
    hpoint = line('XData',x(1),'YData',y(1),'ZData',z(1),'Color','black','Marker',...
        'o','MarkerSize',10,'MarkerFaceColor','black');
    ht = hgtransform('parent',ah);
    set(hpoint,'Parent',ht);

    for i=2:length(x)
        tx = x(i)-x(i-1);
        ty = y(i)-y(i-1);
        tz = z(i)-z(i-1);
        trans = makehgtform('translate',[tx ty tz]),      
        set(ht,'Matrix',trans);
        pause(0.01);
    end