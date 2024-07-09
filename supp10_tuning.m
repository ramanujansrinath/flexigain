clc; clear; close all;
load('data/supp10_tuning_shapes.mat');
load('data/supp10_tuning_resp.mat')

%%
colorDiff = (mod([params.stimRF_num],5)-mod([params.stimOpp_num],5)) == 0;
eid = [3*ones(1,32) ones(1,64) 3*ones(1,32)];
resp = resp(:,eid==1);
resp_base = resp_base(:,eid==1);
good = mean(resp)>1.1*mean(resp_base); % true(1,size(resp,2));
resp = resp(:,good);

figure('color',[0.2 0.2 0.2],'pos',[171,359,356,327]);
for cc=1:47
    clf;
    rr = groupsummary(resp(:,cc),[params.stimRF_num]','mean');
    rr = (rr-min(rr))/(max(rr)-min(rr));

    ha = tight_subplot(nCol,nShape,0.05);
    ha = reshape(ha,nShape,nCol)';
    for jj=1:nCol
        for ii=1:nShape
            % subplot(nCol,nShape,ii+(jj-1)*nShape)
            patch(sh{ii}(:,1),sh{ii}(:,2),cols(jj,:),'facealpha',1,'edgecolor','none','parent',ha(jj,ii));
            set(ha(jj,ii),'color',[rr(5*(jj-1)+ii) 0 0])
            % axis(ha(jj,ii),[-1.1 1.1 -1.1 1.1])
            axis(ha(jj,ii),'equal'); % axis(ha(jj,ii),'off');
        end
    end
    pause;
end

