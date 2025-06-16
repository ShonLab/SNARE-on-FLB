function [nmol, xpos, ypos, A, meanI] = particle_detection(img, threshold, dispres)

    % Preprocess image
    img = double(img); % Convert to double
    img_filtered=img-medfilt2(img,[21, 21]); 
    img_filtered=imfilter(img_filtered,fspecial("average",3));

    % particle detection using intensity threshhold
    img_filtered(img_filtered<threshold)=0;
    img_filtered(1:10,:)=0; img_filtered(end-9:end,:)=0; img_filtered(:,1:10)=0; img_filtered(:,end-9:end)=0; 
    final_img=img_filtered;    
    labeled_img=bwlabel(final_img);

    %calculate properties of detected particles (x,y,area,intensity)
    stats = regionprops(labeled_img, 'PixelIdxList');

    xpos=[]; ypos=[]; A=[]; meanI=[];
    for n=1:numel (stats)
        % Extract pixels of the current label
        pixels = stats(n).PixelIdxList;
        if numel(pixels)<5 continue; end

        % Create a sub-image for the current label
        sub_img = false(size(img_filtered));
        sub_img(pixels) = true;

        % Fill the holes in the sub-image
        sub_img_filled = imfill(sub_img, 'holes'); % Fill empty regions

        % Get newly filled pixels (difference between filled and original)
        new_pixels = find(sub_img_filled & ~sub_img);

        % Combine original and newly filled pixels
        all_pixels = union(pixels, new_pixels); % Merge the indices

        % Create an updated sub-image with all pixels
        sub_img_combined = false(size(img_filtered));
        sub_img_combined(all_pixels) = true;


        % Find regional maxima within the combined sub-image
        regional_max = imregionalmax(imgaussfilt(double(sub_img_combined.*img), 2), 8);

        % Get coordinates of regional maxima
        [max_y, max_x] = find(regional_max);
      
        if ~isempty(max_x)
            [y,x]=ind2sub(size(img),all_pixels);
            if numel(max_x)>numel(x)
                continue;
            end
            %for seperation if multiple maxima exist
            cluster=kmeans([y,x],numel(max_x),'start',[max_y,max_x]);
            for k=1:numel(max_x)
                yk=y(cluster==k); xk=x(cluster==k);
                A=[A;numel(xk)];
                pixels=sub2ind(size(img),yk,xk);  
                meanI=[meanI; mean(img(pixels),'all')];

                [~, max_idx] = max(img(sub2ind(size(img), yk, xk)));
                max_yk = yk(max_idx);
                max_xk = xk(max_idx);
                cnt=centroid(img(max_yk+(-2:2),max_xk+(-2:2)));
                
                xpos = [xpos; max_xk+cnt(2)-3];
                ypos = [ypos; max_yk+cnt(1)-3];
            end
        end
    end
    nmol=numel(xpos);

    % 6. Plot results
    if dispres
        figure;
        ax1 = subplot(131);
        imshow3(img, 10);
        title('raw Image');

        ax1 = subplot(132);
        imshow3(img_filtered, 10);
        title('Filtered Image');
        
        ax2 = subplot(133);
        imshow3(img, 10); hold on;
        plot(xpos, ypos, 'yo', 'markersize', 10);
        title(['Detected molecules: ', num2str(numel(xpos))]);
    end
end
