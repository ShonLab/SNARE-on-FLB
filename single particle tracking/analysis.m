clear all; close all;
%% adjustable
% experimental condition
framerate=0.03;  %second
scale=216*10^(-9);  %meter
sz=1; %diameter of particle, unit:pixel 

dispres=0; 
threshold=10; %orgin: 10
target_frame=15;

r0=200; c0=180;
nr=200; nc=200;

len=5; %minimum length that the particle reside in ROI, origin: 5
maxdisp=7; %the maximum distance particles can move after one frame

path={'raw'};
labels={'4 nM','8 nM','15 nM','30 nM'};

for p=1:size(path,2)
    datainfo=dir(path{p});
    nfile{p}=length(datainfo)-2;
    r_folder=strcat('result_',path{p}); %folder to save results
    mkdir(r_folder)
    for i=1:nfile{p}
            D_data{p}{i}=[]; alpha_data{p}{i}=[];
            spot_I{p}{i}=[];

           fname{p}{i}=datainfo(i+2).name;
           if contains(fname{p}{i},'0. ') continue; end
            if contains(fname{p}{i},'%1600') continue; end
           img=strcat(path{p},'\',fname{p}{i});
           info=imfinfo(img);
           for imageNumber=1:size(info,1)
               img_raw(:,:,imageNumber)=double(imread(img,'index',imageNumber));
           end
           bkg=medfilt2(min(img_raw,[],3),[30, 30]);
           filtered_img=img_raw(r0:r0+nr-1,c0:c0+nc-1,:)-bkg(r0:r0+nr-1,c0:c0+nc-1);
        
           tr{p}{i}=particle_tracking(filtered_img,threshold, maxdisp, len);
           if isempty(tr{p}{i}) continue; end
           rname=strcat(r_folder,'\',fname{p}{i}); delete(rname);
           traj=trajectory(size(img_raw(r0:r0+nr-1,c0:c0+nc-1,:)),tr{p}{i});  
           writeTIFrgb(img_raw(r0:r0+nr-1,c0:c0+nc-1,:),traj,zeros(size(filtered_img)),[r_folder,'\traj_',fname{p}{i}]);

           particle_list=unique(tr{p}{i}(tr{p}{i}(:,3)==target_frame,4));
           D_list=CalD(2,tr{p}{i}(ismember(tr{p}{i}(:,4),particle_list),:), framerate, 4,scale,rname);
           if isempty(D_list)  continue; end  

           tr_roi=tr{p}{i}(tr{p}{i}(:,3)==target_frame,:);
           tr_roi = tr_roi(ismember(tr_roi(:,4), D_list(:,3)), :);
           for mol=1:size(tr_roi,1)
                r=round(tr_roi(mol,2)+(-sz:sz)); c=round(tr_roi(mol,1)+(-sz:sz));
                spot_I{p}{i}(mol,1)=mean(filtered_img(r,c,target_frame),'all');
                D_data{p}{i}(mol,1)=D_list(find(D_list(:,3)==tr_roi(mol,4)),1)*10^12;
                alpha_data{p}{i}(mol,1)=D_list(find(D_list(:,3)==tr_roi(mol,4)),2);
           end

        clear filtered_img img_raw traj
    end

end

save('data.mat')

%% data loading
clear all; close all;
load('data.mat')

dat1={}; dat2={}; dat3={};
count=1;
temp1=[]; temp2=[]; temp3=[]; temp_g=[];
for idx=50:-1:11
    temp1=vertcat(temp1,cell2mat(D_data{1}(idx)));
    temp2=vertcat(temp2,cell2mat(alpha_data{1}(idx)));
    temp3=vertcat(temp3,cell2mat(spot_I{1}(idx)));

    if rem(idx,10)==1
    idx=1:numel(temp2);
    dat1{count}=temp1; dat2{count}=temp2; dat3{count}=temp3;
    

    temp1=[]; temp2=[]; temp3=[];
    count=count+1;
    end
end
%% diffusion coefficient and exponent

figure;
offset=0.2;
yyaxis left
for e=1:4
    plotSpread(dat1{e},'xValues',e-offset, 'SpreadWidth',0.4, 'distributionColors','b')
    boxplot(dat1{e},'Positions',e-offset,'PlotStyle', 'compact','Colors','k','Symbol','');
end
ylabel('log(D [μm^2/s])')

yyaxis right
for e=1:4
    plotSpread(dat2{e},'xValues',e+offset, 'SpreadWidth',0.4, 'distributionColors','r')
    boxplot(dat2{e},'Positions',e+offset,'PlotStyle', 'compact','Colors','k', 'Symbol','');
end
ylabel('diffusion exponent')
xticks(1:4)
xticklabels(labels)
xlabel('concentration (nM)')

%% spot intensity
figure;
g=[];
for e=1:4
    g=cat(1,g,cellstr(repmat(labels{e},[numel(dat3{e}),1])));
end

plotSpread(dat3,'SpreadWidth',0.4,'distributionColors',[0.5 0.5 0.5]); hold on
h=boxplot(cell2mat(dat3'),g,'Colors','k', 'Symbol',''); hold off
set(h, 'LineWidth', 2);
ylim([15,100])
ylabel('spot intensity of a647 (a.u.)')
xticks(1:4)
xticklabels(labels)
xlabel('concentration (nM)')