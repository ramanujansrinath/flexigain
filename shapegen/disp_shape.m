function disp_shape(shape,showAxis)
    h = gca; cla(h);
    if ~exist('showAxis','var'); showAxis = 0; end

    vert = shape.vert;
    face = shape.face;
    
    hp = patch('vertices',vert,'faces',face,'parent',h); hold(h,'on');
    hp.EdgeColor = 'none'; 
    hp.FaceColor = shape.params.color;
    hp.FaceAlpha = 1;
    hp.FaceLighting = 'gouraud';
    material(hp,[0.2 0.7 shape.params.gloss 20 0.5]);
    
    if showAxis
        ax = shape.axis;
        hp.FaceAlpha = 0.7;
        plot3(ax(:,1),ax(:,2),ax(:,3),'-','LineWidth',3,'color','r'); % 1-shape.params.color
    end

    ch = get(h,'children');
    lightExists = sum(arrayfun(@(x) contains(class(ch(x)),'Light'),1:length(ch)));
    if ~lightExists
        light('parent',h,'Position',shape.params.light); 
    end
    % shape.params.light(3) = -shape.params.light(3);
    % light('parent',h,'Position',shape.params.light);
    % hl.Style = 'infinite';
    
    axis(h,'off'); axis(h,'equal');
    set(h,'Projection','perspective')
    set(h,'CameraUpVector',[0 1 0])
    set(h,'CameraViewAngle',40)
    set(h,'CameraPosition',[0,1,5]);
    set(h,'CameraTarget',mean(shape.vert));
end