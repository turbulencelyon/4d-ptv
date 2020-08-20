function [X,Y,Z]=visu_1(tracks,elm)

        X=tracks(elm).x;
        Y=tracks(elm).y;
        Z=tracks(elm).z;
        plot3(X,Y,Z);
        hold on
xlabel('x');
ylabel('z');
zlabel('y');

end
