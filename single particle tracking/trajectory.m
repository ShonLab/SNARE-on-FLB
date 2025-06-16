function out=trajectory(sz, tr)
%function for calculation trajectories of particles
%input:
%   tr: result of the function "track"
%       (x y t id) matrix, which is sorted according to the time
%output:
%   out: trajectory img

ppre=[];    %x y t id
img=zeros([sz(1) sz(2)]);
    for imageNumber = 1:sz(3)
        [r,~]=find(tr(:,3)==imageNumber);  %particle rows for partices that exists in the frame at time t
        line=[];
        for p=1:length(r)   %particle
            if isempty(ppre) || isempty(find(ppre(:,4)==tr(r(p),4)))    
                ppre=vertcat(ppre, tr(r(p),:));
            else 
                [rp, ~]=find(ppre(:,4)==tr(r(p),4));   
                [x,y]=bresenham(ppre(rp, 1),ppre(rp, 2), tr(r(p),1), tr(r(p),2));
                img(sub2ind(size(img),y,x))=1;
                ppre(rp,:)=tr(r(p),:);
            end
        end
        out(:,:,imageNumber)=img;
    end
end

