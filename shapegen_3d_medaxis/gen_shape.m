function shape = gen_shape()
    % medial axis
    params.ma.pts = 31;
    params.ma.curve = rand;
    params.ma.length = 3.5+1*rand;
    params.ma.axis = getCurvedAxis(params.ma);

    % surface
    params.surface.cs.comp = rand;
    [cPts,dPts] = gen_shape_2d(params.surface.cs.comp,false);
    params.surface.cs.cPts = cPts;
    params.surface.cs.dPts = dPts;
    [radProf,radProfXY] = getRadProf(params.ma.pts,0.1+0.6*rand(1,3),[]);
    params.surface.radProfXY = radProfXY;
    params.surface.radProf = radProf;
    params.surface.twistProf = linspace(0,rand*pi,params.ma.pts); 
    params.surface.doSmooth = 1;

    % global
    params.pos = [rand(1,2) 0];
    params.rot = [pi.*rand (pi/2-pi/6)+(pi/3*rand) pi.*rand];
    params.size = 0.4+0.6*rand;
    params.color = getColor;
    params.gloss = 0.001+0.7*rand;
    tt = 0; % -10 + 20*rand;
    params.light = [tt; 5; sqrt(100-(tt.^2))];
    
    % final gen shape
    params.morphType = 0;
    [vert,face,ax] = makeShape(params);
    shape.vert = vert;
    shape.face = face;
    shape.axis = ax; % this is because the axis in shape.params.ma is not the right size or pos
    shape.params = params;
end

function [vert,face,ax] = makeShape(params)
    cs = params.surface.cs.dPts;
    nCsPts = size(cs,1);
    nVert = size(cs,1)*params.ma.pts - nCsPts+1 - nCsPts+1;
    
    vert = nan(nVert,3);
    
    % vert - assemble cs
    for ii=1:params.ma.pts
        csspin = movePts(cs,0,0,1,rad2deg(params.surface.twistProf(ii)),[0 0]);
        cssc = [csspin(:,1).*params.surface.radProf(ii,1) csspin(:,2).*params.surface.radProf(ii,2) zeros(nCsPts,1)]; % scale
        csrot = rotvert(cssc,[params.ma.axis.ang(ii) 0 0]); % rotate for ma
        
        % this is hacky: just a way for the capping points to be closer
        % to shape to make it look smoother.
        % if ii==1
        %     pos = (params.ma.axis.pos(ii+2,:)); % +params.ma.axis.pos(ii,:)) ./ 2;
        % elseif ii==2
        %     pos = (params.ma.axis.pos(ii+1,:)); % +params.ma.axis.pos(ii,:)) ./ 2;
        % elseif ii==params.ma.pts-1
        %     pos = (params.ma.axis.pos(ii-1,:)); % +params.ma.axis.pos(ii,:)) ./ 2;
        % elseif ii==params.ma.pts
        %     pos = (params.ma.axis.pos(ii-2,:)); % +params.ma.axis.pos(ii,:)) ./ 2;
        % else
            pos = params.ma.axis.pos(ii,:);
        % end
        
        sec =  csrot + repmat(pos,nCsPts,1); % shift

        if ii==1
            vert(1,:) = sec(1,:);
        elseif ii==params.ma.pts
            vert(nVert,:) = sec(1,:);
        else
            vert(nCsPts*(ii-2)+2 : 1+nCsPts*(ii-1), :) = sec;
        end

        % plot3(vert(:,1),vert(:,2),vert(:,3),'r.'); 
        % axis equal; axis off; set(gca,'XLim',[-2 2]); view(-14,-14);
    end
    ax = params.ma.axis.pos;
    ax = ax - repmat((min(vert)+max(vert))/2,size(ax,1),1);
    vert = vert - repmat((min(vert)+max(vert))/2,nVert,1);
    vert = round(vert,4);

    filepath = fileparts(mfilename('fullpath'));
    facepath = [filepath '/export/face_' num2str(params.ma.pts) '.txt'];
    if ~exist(facepath,'file')
        face = [];
        % face - end cap
        for ii=2:nCsPts+1
            face = [face ; 1 ii (ii+1)]; %#ok<AGROW>
        end
        face(end,3) = 2;

        % face - middle
        for ii=2:params.ma.pts-2
            thisRow = 2+(ii-2)*nCsPts : 1+(ii-1)*nCsPts;
            nextRow = 2+(ii-1)*nCsPts : 1+ii*nCsPts;
            % plot3(vert(thisRow,1),vert(thisRow,2),vert(thisRow,3),'r.'); hold on
            % plot3(vert(nextRow,1),vert(nextRow,2),vert(nextRow,3),'b.')

            for jj=1:nCsPts
                if jj<nCsPts
                    face = [face; thisRow(jj) nextRow(jj) nextRow(jj+1)]; %#ok<AGROW>
                    % patch(vert(face(end,:),1),vert(face(end,:),2),vert(face(end,:),3),'c','facealpha',0.2,'edgecolor','none');
                    face = [face; thisRow(jj) nextRow(jj+1) thisRow(jj+1)]; %#ok<AGROW>
                    % patch(vert(face(end,:),1),vert(face(end,:),2),vert(face(end,:),3),'c','facealpha',0.2,'edgecolor','none');
                else
                    face = [face; thisRow(jj) nextRow(jj) nextRow(1)]; %#ok<AGROW>
                    % patch(vert(face(end,:),1),vert(face(end,:),2),vert(face(end,:),3),'k','facealpha',0.8);
                    face = [face; thisRow(jj) nextRow(1) thisRow(1)]; %#ok<AGROW>
                    % patch(vert(face(end,:),1),vert(face(end,:),2),vert(face(end,:),3),'k','facealpha',0.8);
                end
            end
        end

        % face - end cap
        for ii=1:nCsPts-1
            face = [face; nextRow(ii) nVert nextRow(ii+1)]; %#ok<AGROW>
        end
        face(end+1,:) = [nextRow(end) nVert nextRow(1)];
        writematrix(face,facepath);
    else
        face = readmatrix(facepath);
    end
    
    ax =  ax .* params.size;
    ax = rotvert(ax,params.rot);
    ax =  ax + repmat(params.pos,size(ax,1),1);
    
    vert =  vert .* params.size;
    vert = rotvert(vert,params.rot);
    vert =  vert + repmat(params.pos,size(vert,1),1);
    
    if params.surface.doSmooth
        vert = smoothVert(vert,face);
    end
end

function vert = smoothVert(vert,face)
    neigh = cell(size(vert,1),1);
    for smoothTimes=1:20
        smVert = zeros(size(vert));
        for ii=1:length(vert)
            if isempty(neigh{ii})
                a = face(sum(face==ii,2)>0,:);
                a = a(:); 
                a(a==ii) = [];
                neigh{ii} = a;
            end
                
            smVert(ii,:) = mean(vert(neigh{ii},:));
        end
        vert = smVert;
    end
end

function vert = rotvert(vert,rot)
    rx = [1 0 0; 0 cos(rot(1)) -sin(rot(1)); 0 sin(rot(1)) cos(rot(1))];
    ry = [cos(rot(2)) 0 sin(rot(2)); 0 1 0; -sin(rot(2)) 0 cos(rot(2))];
    rz = [cos(rot(3)) -sin(rot(3)) 0; sin(rot(3)) cos(rot(3)) 0; 0 0 1];

    vert = vert * rx;
    vert = vert * ry;
    vert = vert * rz;
end

function radProfXY = getRadProfXY(surfThkns)
    if size(surfThkns,1) == 5
        radProfXY = surfThkns;
    else
        x = [0 0.1+0.2*rand 0.4+0.2*rand 0.7+0.2*rand 1];
        y = [0 surfThkns 0];
        radProfXY = [x' y'];
    end
end

function [rProf,cPts] = getRadProf(nPts,surfThkns,surfProf)
    if ~isempty(surfProf)
        rProf = surfProf;
        cPts = [];
        return
    end

    tryTimes = 0;
    maxTryTimes = 100;
    
    % OLD
    % while tryTimes < maxTryTimes
    %     cPts = getRadProfXY(surfThkns);
    %     cPts = [cPts(1,:); 0 mean(cPts(1:2,2)); cPts(2:end,:)];
    %     cPts(end+1:end+2,:) = [cPts(end,:);cPts(end,:)];
    %     cPts(end-2,2) = mean(cPts(end-3:end-2,2));
    % 
    %     dPts = getDensePoints(cPts,[],nPts,0);
    %     yy = interp1(dPts(:,1),dPts(:,2),linspace(0,1,nPts),'spline');
    %     yy = yy - max(yy([1 nPts]));
    %     yy(1) = 0;
    %     yy(nPts) = 0;
    %     % yy = interp1(linspace(0,1,nPts),yy,linspace(0,1,nPts),'spline');
    %     if sum(yy<0)
    %         tryTimes = tryTimes + 1;
    %     else
    %         break;
    %     end
    % end
    
    % simpler but pointy ended
    while tryTimes < maxTryTimes
        cPts = getRadProfXY(surfThkns);
        xx = linspace(0,1,nPts);
        yy = interp1(cPts(:,1),cPts(:,2),xx,'spline');
        if sum(yy<0)
            tryTimes = tryTimes + 1;
        else
            break;
        end
    end
    
    if tryTimes == maxTryTimes
        error('Couldn''t find a good rad profile');
    end
    
    rProf(:,1) = yy;
    rProf(:,2) = rProf(:,1);
    rProf = round(rProf,4);
end
    
function ax = getCurvedAxis(ma)
    len = ma.length;
    curve = ma.curve;
    nPts = ma.pts;
    
    if ma.curve == 0
        rad = 1/0.000001;
    else
        rad = 1/curve; 
    end
    
    % uneven sampling
    % len = len/pi; % specified len is arc length, not rad
    % x = linspace(-len/2,len/2,nPts);
    % y = sqrt(rad^2-x.^2);
    % y = y-max(y);
    
    th = (pi+len/rad)/2 - linspace(0,len/rad,nPts);
    th(end+1) = th(end)-(th(end-1)-th(end)); % for tangents
    x = rad*cos(th);
    y = rad*sin(th);
    y = y-max(y);

    ma = [zeros(1,nPts+1); y; x]';
    tang = circshift(ma,-1) - ma;
    ma(end,:) = [];
    tang(end,:) = [];
    tang = tang./repmat(sqrt(sum(tang.^2,2)),1,3); % unit vec
    
    ax.pos = ma;
    ax.tan = tang;
    ax.ang = atan2(tang(:,2),tang(:,3));
    
    % plot3(ma(:,1),ma(:,2),ma(:,3),'b.'); axis equal; grid on; view(0,90); hold on;
    % for ii=1:nPts
    %     p1 = ma(ii,:);
    %     p2 = p1 + 0.2*tang(ii,:);
    %     hl = line([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'color','r','linewidth',2);
    %     pause
    %     delete(hl)
    % end
end
