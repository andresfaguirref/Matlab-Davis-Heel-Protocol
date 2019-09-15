function Iter = Iter_point(point,frames,plot)

N_point = [];
N_frames = [];
y=1;
for x=1:length(point)
    if not(isnan(point(x)))
        N_point(y) = point(x);
        N_frames(y) = frames(x);
        y = y+1;
    end
end

Iter = spline(N_frames,N_point,frames);
if plot 
    figure()
    plot(N_frames,N_point,'o',frames,Iter,'*r')
end

end