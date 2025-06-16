function [img_1,img_2]=dualviewer_spliter(movie,tform)
    nframe=size(movie,numel(size(movie)));
    for imageNumber=1:nframe
        img_2(:,:,imageNumber)=movie(1+end/2:end,:,imageNumber);

        img_1_temp = imwarp(movie(1:end/2,:,imageNumber), tform, 'OutputView', imref2d(size(img_2(:,:,imageNumber))));
        img_1(:,:,imageNumber)=img_1_temp;     
    end
end
