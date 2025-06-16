clear all; close all;
%% adjustable
maxdisp=7; %the maximum distance particles can move after one frame

r0=22; c0=110;
nr=200 ; nc=200 ;

%parameters for colocalization
target_frame=25;

thres=[15 15];
disp=[100, 150]; 

max_off=4; 
%%
plist=dir;
isSubfolder=[plist.isdir]; isSubfolder(1:2)=0;
path={plist(isSubfolder).name};
path(contains(path, 'result'))=[];
path(contains(path, 'code'))=[];
load('tform.mat'); 

for p=1:size(path,2)
% for p=2
    r_folder=['result_' path{p}]; %folder to save results
    mkdir(r_folder)

    datainfo=dir(path{p});
    fname{p}={datainfo.name};
    fname{p}(~contains(fname{p},'.tif'))=[];

    %% bkg calculation
    count=1;
    bkg_idx=find(contains(fname{p},'background'));
    for i=bkg_idx
        img=strcat(path{p},'\',fname{p}{i});
        info=imfinfo(img);
        for imageNumber=1:size(info,1)
            img_raw(:,:,imageNumber)=imread(img,'index',imageNumber);
        end
        % bkgs for individual movies 
        bkg_temp(:,:,count)=min(img_raw,[],3);
        count=count+1;
    end
    bkg=uint16(medfilt2(mean(bkg_temp,3),[20,20]));

    fname{p}(bkg_idx)=[];
    clear img_raw bkg_temp

    %%
    exps{p}=cellfun(@(x) x(1), fname{p});
    exp_labels=unique(exps{p});
    for e=1:numel(exp_labels)
        nexp{p}(e)=sum(exps{p}==exp_labels(e));
    end
    nfile{p}=numel(fname{p});

    for i=1:nfile{p}
       colocal{p}{i}=[]; 
       spot_I1{p}{i}=[]; 

       img=strcat(path{p},'\',fname{p}{i});
       info=imfinfo(img);

       %image loading
       for imageNumber=1:30
%        for imageNumber=1:5
           img_raw(:,:,imageNumber)=imread(img,'index',imageNumber);
       end

       if i==find((exps{p}==fname{p}{i}(1)),1)
           background=bkg;
           leakage=0;
%            [tform,comp_img]=find_tform(img_raw(:,:,1)-background,disp);
       end

       %colocalization
       [up,down]=dualviewer_spliter(img_raw-background,tform);   
       down=down-up*leakage; %leakage substraction
       [co_pos{p}{i},~]=colocalization(strcat(r_folder,'\colocal_',fname{p}{i},'.png'),up(r0:r0+nr-1,c0:c0+nc-1,target_frame), down(r0:r0+nr-1,c0:c0+nc-1,target_frame),thres,max_off,disp);

       writeTIFrgb(down,zeros(size(up)),up,strcat(r_folder,'\composite_',fname{p}{i}));
       
       %% spot intensity analysis
       [nmol{p}{i},xpos,ypos,~,~]=particle_detection(up(r0:r0+nr-1,c0:c0+nc-1,target_frame),thres(1),0);

       for mol=1:nmol{p}{i}
           r=round(ypos(mol)+(-2:2)); c=round(xpos(mol)+(-2:2));              
           spot_I1{p}{i}(mol,:)=[mean(up(r0+r-1,c0+c-1,target_frame),'all'),mean(down(r0+r-1,c0+c-1,target_frame),'all')];
       end

        %% spot intensity analysis_DID
       [nmol_temp,xpos,ypos,~,~]=particle_detection(down(r0:r0+nr-1,c0:c0+nc-1,target_frame),thres(2),0);
       mean_I{p}{i}=mean(min(down(50:150,150:250,:),[],3),"all");

       for mol=1:nmol_temp
           r=round(ypos(mol)+(-2:2)); c=round(xpos(mol)+(-2:2));              
           spot_I2{p}{i}(mol,:)=[mean(up(r0+r-1,c0+c-1,target_frame),'all'),mean(down(r0+r-1,c0+c-1,target_frame),'all')];
          if isempty(cell2mat(co_pos{p}{i}{2})) colocal{p}{i}(mol)=0; 
           else colocal{p}{i}(mol)=sum(ismember([xpos(mol) ypos(mol)],cell2mat(co_pos{p}{i}{2}),'rows'));
           end
       end
       %%
       clear img_raw up down mov
    end
end

save(['data.mat'])


%% data loading
clear
load data.mat
%%
figure;
for p=1:numel(path)
    temp{p}=cell2mat(spot_I1{p}');
end

figure;
edges=0:50:1800;
histogram(temp{1}(:,1),edges,'Normalization','probability','EdgeAlpha',0); hold on;
histogram(temp{2}(:,1),edges,'Normalization','probability','EdgeAlpha',0); hold off;
legend({'wo HD','w HD'})
ylabel('probability')
xlabel('Egfp intensity (a.u.)')
%%
figure;
for p=1:numel(path)
    temp{p}=cell2mat(spot_I2{p}');
end

figure;
edges=0:50:1800;
histogram(temp{1}(:,1),edges,'Normalization','probability','EdgeAlpha',0); hold on;
histogram(temp{2}(:,1),edges,'Normalization','probability','EdgeAlpha',0); hold off;
legend({'wo HD','w HD'})
ylabel('probability')
xlabel('DiD intensity (a.u.)')

%%
for p=1:numel(path)
    temp{p}=[];
    for i=1:nfile{p}
        temp{p}=[temp{p}; mean_I{p}{i}];
    end
end

figure(Position=[500 300 200 400]);
boxplot(temp{1},'Positions',1); hold on
boxplot(temp{2},'Positions',2); hold on

plotSpread({temp{1},temp{2}},'distributionColors',{'b','r'},'xNames',{'wo HD','w HD'}); 
ylabel('mean background intensity')
hold off
%%
figure(Position=[500 300 200 400]);
for p=1:numel(path)
    temp{p}=[];
    for i=1:nfile{p}
        temp{p}=[temp{p}; sum(colocal{p}{i},"all")/numel(colocal{p}{i})];
    end
end
boxplot(temp{1},'Positions',1); hold on
boxplot(temp{2},'Positions',2); hold on

plotSpread({temp{1},temp{2}},'distributionColors',{'r','b'},'xNames',{'w HD','wo HD'}); 
ylabel('DiD to Syn ratio')
