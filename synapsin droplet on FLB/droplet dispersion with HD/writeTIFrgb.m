function writeTIFrgb(R,G,B,fname)
    delete(fname)

    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.Photometric = Tiff.Photometric.RGB;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 3;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    
    if numel(size(R))==3
        nframe=size(R,3);
    else
        nframe=1;
    end

    
    for i=1:nframe
        tagstruct.ImageLength = size(R(:,:,i), 1);
        tagstruct.ImageWidth = size(R(:,:,i), 2);
        composite=Tiff(fname,'a');
        composite.setTag(tagstruct);
        comp_img=cat(3,R(:,:,i),G(:,:,i),B(:,:,i));
        composite.write(uint16(comp_img));
    end
       
    close(composite)
end
