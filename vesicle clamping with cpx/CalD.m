function out=CalD(dim, tr, framerate,len, scale,fname)
% calculate the diffusion coefficient using the formula:
% msd=2nDt-->D=msd/2nt
%input:
%   dim: the dimension of the area where particles is moving in
%   lineM: line matrix, the result of the function "trajectory"
%          (x; y; t; id)
%   framerate: in unit 'sec'
%   scale: the pixel size, in unit 'meter'
%
%output:
%   out: the list of the estimated diffusion coefficients of the each particle
%        (D id)
%        in unit m^2/s


    tr=sortrows(tr,4);    %sort according to its id
    particles=unique(tr(:,end));  %particle list
    msd={}; tau={};
    D=[];
    error_thres=0.1;
    
    figure;
    for p=1:size(particles,1)
        [r,~]=find(tr(:,4)==particles(p));
        if isempty(r) break; end
        N=size(r,1);
        x=tr(r,1); y=tr(r,2);
        for dt=2:N
            dx=x(dt:end)-x(1:end-dt+1);
            dy=y(dt:end)-y(1:end-dt+1);
            msd{p}(dt-1)=mean(dx.^2+dy.^2)*scale^2;
            tau{p}(dt-1)=(dt-1)*framerate;
        end
        
        log_msd=log10(msd{p}(1:len)); log_tau=log10(tau{p}(1:len));
        coeffs=polyfit(log_tau,log_msd,1);
        
        R2=1-(sum((log_msd-polyval(coeffs, log_tau)).^2)/sum((log_msd-mean(log_msd)).^2));

        if R2>1-error_thres
            loglog(tau{p},msd{p}); hold on;
            D=[D; 10^(coeffs(2))/2/dim, coeffs(1), particles(p)];
        else
            loglog(tau{p},msd{p},'k'); hold on;
        end   
    end

    hold off;

    ylabel('msd (μm^2)')
    xlabel('tau (s)')
    saveas(gcf,[fname,'.fig'])
    
    out=D;
end

