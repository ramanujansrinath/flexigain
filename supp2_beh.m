clc; close all; clear;
load('data/supp2_beh.mat','beh_main')
setList = 1:120;

%%
figure('color','w','position',[53,216,2508,787])
ha = tight_subplot(7,21,[0.02 0.01],[0.04 0.01],[0.03 0.01]); 
ha = reshape(ha,21,7)';

imgPath = 'data/supp2_img';

% show stimulus images on top
for jj=1:20
    % load([limbPath '/img_' num2str(setList(jj)) '_5.mat'])
    % imshow(theImage,'parent',ha(1,1+jj)); 
    imshow([imgPath '/img_' num2str(setList(jj)) '_5.png'],'parent',ha(1,1+jj)); 
end

% show orientation and color changes on the left
yim = [1 21 41 61 81 101]; anno = {'0{\circ}' '45{\circ}' '90{\circ}' '135{\circ}' 'alt color' 'gray'};
for jj=1:length(yim)
    % load([limbPath '/img_' num2str(setList(jj)) '_5.mat'])
    % imshow(theImage,'parent',ha(1,1+jj)); 
    imshow([imgPath '/img_' num2str(yim(jj)) '_5.png'],'parent',ha(1+jj,1)); 
    ht = text(1,1,anno{jj},'fontname','lato','parent',ha(1+jj,1),'Interpreter','tex');
end
delete(ha(1,1));

% plot behavior per monkey
monkeys = {'ivory' 'pattern'};
cols = [0.9 0.5 0.2; 0.2 0.5 0.9];
for mm=1:2
    beh_monk = beh_main(cellfun(@(x) contains(x,monkeys{mm}),{beh_main.monkey}));
    
    for jj=1:length(setList)
        if mod(setList(jj),20)==0
            rowNum = 1+floor(setList(jj)/20);
            colNum = 21;
        else
            rowNum = 2+floor(setList(jj)/20);
            colNum = 1+mod(setList(jj),20);
        end
        h = ha(rowNum,colNum); hold(h,'on');
        binEdges = linspace(0,1,9);
        curv_m = (binEdges+circshift(binEdges,-1))/2; curv_m = curv_m(1:end-1);

        curv = [beh_monk([beh_monk.set] == setList(jj)).curv];
        if ~isempty(curv)
            sel = [beh_monk([beh_monk.set] == setList(jj)).sel];
            [~,~,groups] = histcounts(curv,binEdges);
            sel_m = groupsummary(sel',groups','mean');
            sel_s = groupsummary(sel',groups','std'); % ./sqrt(groupsummary(sel',groups','nnz'));
            
            patch([curv_m'; flipud(curv_m')],[(sel_m - sel_s/2);(flipud(sel_m) + flipud(sel_s)/2)],cols(mm,:),'edgecolor','none','facealpha',0.2,'parent',h)
            % errorbar(h,x,y,s,']color','k','capsize',0,'linewidth',2)
            plot(h,curv_m,sel_m,'-','markersize',6,'markerfacecolor','w','linewidth',2,'color',cols(mm,:))
        end

        if mm==2
            line([0 1],[0 1],'color',[0.7 0.7 0.7],'linestyle','--','linewidth',2,'parent',h)
            % text(0,1,num2str(setList(jj)),'FontSize',15,'FontName','lato','Color',[0.5 0.5 0.5],'Parent',h)
            set(h,'XTickLabel',0:0.5:1,'YTickLabel',0:0.5:1)
            fixPlot(h,[-0.1 1.1],[-0.1 1.1],'','',0:0.5:1,0:0.5:1);
            grid(h,'off');
    
            if mod(jj,20)==1 && jj<101
                set(h,'xtick',[])
            elseif jj>101
                set(h,'ytick',[])
            elseif jj<101
                set(h,'xtick',[],'ytick',[])
            end
        end
    end
end
