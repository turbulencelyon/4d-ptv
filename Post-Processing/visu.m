function visu(tracks,seuil,color)



index=find([tracks.L]>seuil);
stop=length(index)
for i=1:stop
elm=index(i);
%test=tracks(elm).nbstitch;
test=1;
    if (~test==0)
        x=tracks(elm).x;
        y=tracks(elm).y;
        z=tracks(elm).z;
        plot3(x,z,y,color);
        hold on

    end 
end
xlabel('x')
ylabel('z')
zlabel('y')

end
