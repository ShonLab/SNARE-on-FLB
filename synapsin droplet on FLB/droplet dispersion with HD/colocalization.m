function [co_pos,co_ratio]=colocalization(fname,img_1,img_2,feature_thres,max_off, disp_range)
    %img_1 for upper, img_2 for below

    if numel(size(img_1))==3 nframe=size(img_1,3);
    else nframe=1; end

    for imageNumber=1:nframe
        if imageNumber==1 draw=1;
        else draw=0; end
%         [nmol1,xpos1,ypos1]=countSM(img_1(:,:,imageNumber),feature_thres(1),draw);
%         [nmol2,xpos2,ypos2]=countSM(img_2(:,:,imageNumber),feature_thres(2),draw);
        [nmol1,xpos1,ypos1,~,~]=particle_detection(img_1(:,:,imageNumber),feature_thres(1),draw);
        [nmol2,xpos2,ypos2,~,~]=particle_detection(img_2(:,:,imageNumber),feature_thres(2),draw);

        co_ratio(imageNumber,1)=0; co_ratio(imageNumber,2)=0;
        if ~nmol1
            co_pos{1}{imageNumber}=[]; co_pos{2}{imageNumber}=[];
            if imageNumber==1; co_ratio(imageNumber,1)=nan; end
        end

        if ~nmol2
            co_pos{1}{imageNumber}=[]; co_pos{2}{imageNumber}=[];
            if imageNumber==1; co_ratio(imageNumber,2)=nan; end
        end

        if ~(nmol1*nmol2)
            continue;
        end

        D=pdist2([xpos1 ypos1],[xpos2 ypos2]);      

        [idx_1,idx_2]=find(D==min(D,[],2) & D==min(D,[],1)& D<=max_off);
        co_pos{1}{imageNumber}=[xpos1(idx_1) ypos1(idx_1)]; 
        co_pos{2}{imageNumber}=[xpos2(idx_2) ypos2(idx_2)]; 
        co_ratio(imageNumber,1)=numel(idx_1)/nmol1; co_ratio(imageNumber,2)=numel(idx_2)/nmol2; 
    end

   %% comp image
    figure;   

   
    img_comp=double(cat(3, img_2(:,:,1)-medfilt2(min(img_2,[],3),round(size(img_2(:,:,1))/10)), img_1(:,:,1)-medfilt2(min(img_1,[],3),round(size(img_1(:,:,1))/10)), zeros(size(img_2(:,:,1))))); 
    img_comp(:,:,1)=img_comp(:,:,1)/disp_range(2); img_comp(:,:,2)=img_comp(:,:,2)/disp_range(1);
    imshow2(img_comp); hold on;
    title(fname);
    if ~isempty(co_pos{2}{1})
        plot(co_pos{2}{1}(:,1),co_pos{2}{1}(:,2),'yo');
    end
    hold off
    saveas(gcf,fname)
end