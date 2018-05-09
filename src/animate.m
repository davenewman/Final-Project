function animate(u,X,Y,tlength)

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
%     [imind,cm] = rgb2ind(im,256);
end