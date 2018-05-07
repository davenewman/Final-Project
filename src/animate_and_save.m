function animate_and_save(u,x,y,tlength)

[X,Y] = meshgrid(x(2:end-1),y(2:end-1));
filename = 'animation.gif';
first_frame = true;
for idx = 1:tlength
    surf(X,Y,u(:,:,idx));
    xlabel('x'),ylabel('y');
    axis([0 2*pi 0 2*pi -100 200]);
    drawnow
    % create gif animation
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if first_frame
      imwrite(imind,cm,filename,'gif','DelayTime',0.04,'Loopcount',inf);
      first_frame = false;
    else
      imwrite(imind,cm,filename,'gif','DelayTime',0.04,'WriteMode','append');
    end
end
