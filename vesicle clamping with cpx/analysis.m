clear all; close all;
%% adjustable
% experimental condition
framerate=0.03;  %second
scale=216*10^(-9);  %meter
sz=1; %diameter of particle 

r0=90; nr=150;
c0=50; nc=150;

dispres=0; 
threshold=[8,50];

len=5;
maxdisp=5; %the maximum distance particles can move after one frame

path1={'0708','0715'};
path2={'raw1','raw2','raw3'};
%%
img='for bkg.tif';
for imageNumber=1:450
   img_raw(:,:,imageNumber)=double(imread(img,'index',imageNumber));
end
bkg=medfilt2(min(img_raw,[],3),[30, 30]);
clear img_raw

%%
for P=1:size(path1,2)
    for p=1:size(path2,2)
        pth=strcat(path1{P},'\',path2{p});
        datainfo=dir(pth);
        nfile{P}{p}=length(datainfo)-2;
        r_folder=strcat(path1{P},'\result_',path2{p}); %folder to save results
        mkdir(r_folder)
    
        %% L532 signal analysis
        for i=1:nfile{P}{p}
            fname{P}{p}{i}=datainfo(i+2).name;          
            if contains(fname{P}{p}{i},'3. ')
               img=strcat(pth,'\',fname{P}{p}{i});
               info=imfinfo(img);
               for imageNumber=1:size(info,1)
                   img_raw(:,:,imageNumber)=double(imread(img,'index',imageNumber));
               end
               
               filtered_img=img_raw-medfilt2(min(img_raw,[],3),[30, 30]); 
               filtered_img2=img_raw-bkg;
               mean_I{P}{p}{i}=mean(filtered_img2(r0:r0+nr-1,c0:c0+nc-1),"all");
               sprintf(img)
               
               tr=particle_tracking(filtered_img(r0:r0+nr-1,c0:c0+nc-1,:),threshold(1), maxdisp, len);
               if isempty(tr) continue; end
               rname=strcat(r_folder,'\',fname{P}{p}{i}); delete(rname);

               delete([r_folder '\traj_' fname{P}{p}{i}])
               traj=trajectory(size(img_raw(r0:r0+nr-1,c0:c0+nc-1,:)),tr);

               writeTIFrgb(img_raw(r0:r0+nr-1,c0:c0+nc-1,:),traj,zeros(size(img_raw(r0:r0+nr-1,c0:c0+nc-1,:))),[r_folder '\traj_' fname{P}{p}{i}]);

               tr_corrected=tr;
               D_list=CalD(2,tr_corrected, framerate,4, scale,rname);
               if isempty(D_list)  continue; end   

                particles=unique(D_list(:,3));

                idx=1;
    
                for particle=1:size(particles,1)
                    rows=find(tr(:,4)==particles(particle));

                    tr_roi=tr(rows,:);
                    target_frame=find(tr_roi(:,3)==10);

                    if isempty(target_frame) continue; end
                    
                    x_range=c0+(round(tr(rows(target_frame),1)-sz:tr(rows(target_frame),1)+sz))-1;
                    y_range=r0+(round(tr(rows(target_frame),2)-sz:tr(rows(target_frame),2)+sz))-1;
                    roi=filtered_img(y_range,x_range,10);
                    roi2=filtered_img2(y_range,x_range,10);

                    I_data{P}{p}{i}(idx,1)=mean(roi,'all');    
                    I_data2{P}{p}{i}(idx,1)=mean(roi2,'all');     
                    D_data{P}{p}{i}(idx,1)=D_list(particle,1)*10^12;
                    alpha_data{P}{p}{i}(idx,1)=D_list(particle,2);
                    idx=idx+1;
                end
            end
            clear filtered_img img_raw tr tr_corrected traj D_list particles roi
        end
    end
end

save('data.mat')
%% plot
clear; 
close all;
load('data.mat')
labels={'Vamp','ΔN+Vamp','ΔN+Vamp+Cpx'};
for p=1:3
    tmp=[];
    for P=1:2
        if p==1 & P==1 continue; end
        tmp=vertcat(tmp,cell2mat([I_data{P}{p}', D_data{P}{p}' alpha_data{P}{p}',I_data2{P}{p}']));
    end
    data{p}=tmp;   
end

for p=1:3
    tmp=[];
    for P=1:2
        if p==1 & P==1 continue; end
        tmp=vertcat(tmp,cell2mat(mean_I{P}{p}'));
    end
    data2{p}=tmp;   
end

%% D
clear mu ft ini
mu=[0.75 1.4 3.7];
eqn=sprintf('a1*exp(-(x-%f)^2/(2*s1^2)) + a2*exp(-(x-%f)^2/(2*s2^2)) + a3*exp(-(x-%f)^2/(2*s3^2))',mu(1),mu(2),mu(3));
ft=fittype(eqn,'independent','x','coefficients',{'a1', 's1', 'a2', 's2','a3', 's3'});

ini{2} = [5, 0.1, 10, 0.25, 3, 0.5]; % vamp on dn

ini{3} = [6, 0.1, 10, 0.25, 10, 0.1]; % vamp, cpx on dn


fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0]);
figure('Position',[300 300 500 400]);

edges=0.1:0.4:30;
x=0:0.01:30;
edge_center=(edges(1:end-1)+edges(2:end))/2;
for p=2:3
    subplot(3,1,p)
    tmpdat=sort(data{p}(:,2));   
    histogram(tmpdat,edges); hold on

    tmpdat(tmpdat>5)=[];

    [N,~]=histcounts(tmpdat,edges); 

    fitResult{p}=coeffvalues(fit(edge_center',N',ft,'StartPoint', ini{p},'lower',zeros(size(ini{p}))));

    sum_graph=[];
    for comp=1:numel(mu)
        sum_graph(comp,:)=fitResult{p}(2*comp-1)*exp(-(x-mu(comp)).^2/(2*fitResult{p}(2*comp)^2));
        plot(x,sum_graph(comp,:)); hold on
    end
    disp(sum(N))
    plot(x,sum(sum_graph),'k'); hold on
    xline(mu,'-r','linewidth',2)
    xlim([0, 8])
    ylim([0,40])
    ylabel(labels{p})
end
xlabel('D (μm^2)')

figure;
for p=2:3
    subplot(3,1,p)
    integral=[];
    for comp=1:numel(mu)
        integral=[integral; sum((fitResult{p}(2*comp-1)*exp(-(x-mu(comp)).^2/(2*fitResult{p}(2*comp)^2))))];
    end
    pie(integral)
%     legend('location','eastoutside')
end

%% fig3-b
figure;
labels={'VAMP','ΔN+VAMP','ΔN+VAMP+Cpx'};
g=[];
for e=1:3
    bar(e,mean(data2{e})); hold on;
    errorbar(e,mean(data2{e}),std(data2{e}))
end

ylabel('mean intensity of Cy3 (a.u.)')
xticks(1:3)
xticklabels(labels)


%% fig3-a
figure('Position',[300 300 500 400]);

g1={}; g2={};
bkg=22;

unit1=23;
unit2=15;


ft= fittype('a * (round(x)>=1)*(poisspdf(round(x), l) / (1 - exp(-l)))','independent', 'x', 'coefficients', {'a','l'});
ini = [80,2];

gaussft=fittype(['a1*exp(-(x-1)^2/(2*s1^2))'], ...
             'independent', 'x', 'coefficients', {'a1','s1'});
gaussini=[20,0.2];


fine=5;
edges=-0.5:1/fine:10.5;
edge_center=(edges(1:end-1)+edges(2:end))/2;
edges2=-0.5:1:10.5;
edge_center2=(edges2(1:end-1)+edges2(2:end))/2;
edges3=-0.5:1/fine/2:10.5;
edge_center3=(edges3(1:end-1)+edges3(2:end))/2;

for p=2:3
    subplot(3,1,p)
    tmpdat=sort(data{p}(:,1))-bkg;
    norm_dat=tmpdat/unit1;

    [N_fine,~]=histcounts(norm_dat,edges); 
    N1=[];
    for i=1:floor(numel(N_fine)/fine)
        N1(i)=mean(N_fine(fine*(i-1)+1:fine*i));
    end

    fitResult{p}=coeffvalues(fit(edge_center2(2:end)',N1(2:end)',ft,'StartPoint', ini));

    g1{p}=fitResult{p}(1)*(round(edge_center2)>=1).*(poisspdf(round(edge_center2), fitResult{p}(2))/(1 - exp(-fitResult{p}(2))));
    g1{p}=interp1(edge_center2,g1{p},edge_center);
    g1{p}(edge_center<1)=0;

    rem=round(sum(N_fine(1:round(1.2*fine))));
    rem_dat=tmpdat(1:rem)/unit2;
    [N2,~]=histcounts(rem_dat,edges); 

    gfitresult=coeffvalues(fit(edge_center(edge_center>0 & edge_center<1.5)',N2(edge_center>0 & edge_center<1.5)',gaussft,'StartPoint', gaussini));
    g2{p}=gfitresult(1)*exp(-(edge_center3-1).^2/(2*gfitresult(2)^2))*unit1/unit2;
    g2{p}=interp1(edge_center3 * unit2 / unit1, g2{p}, edge_center);

    histogram(norm_dat*unit1, edges*unit1); hold on    
    
    plot(fine_x * unit1, g1_interp, 'y'); hold on 
    plot(fine_x * unit1, g2_interp, 'r'); hold on;
    plot(fine_x * unit1, g1_interp+g2_interp, 'k'); hold on;

    xlim([0,150])

    ylabel(labels{p})
    ylim([0,40])
end
xlabel('spot intensity of Cy3')

figure;
for p=2:3
    subplot(3,1,p)
    pie([sum(g2{p}(~isnan(g2{p}))); sum(g1{p}(~isnan(g1{p})))])
    legend('location','eastoutside')
end
