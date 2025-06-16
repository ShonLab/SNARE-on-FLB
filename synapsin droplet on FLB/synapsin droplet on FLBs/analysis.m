clear all; close all;
%% adjustable
r0=15; c0=115;
nr=220 ; nc=220 ;

%parameters for colocalization
target_frame=1;

thres=[15 20];
disp=[100, 500]; 

max_off=4; 

%%
plist=dir;
isSubfolder=[plist.isdir]; isSubfolder(1:2)=0;
path={plist(isSubfolder).name};
path(contains(path, 'result'))=[];
load('tform.mat'); load('bkg.mat')
for p=1:size(path,2)
    r_folder=['result_' path{p}]; %folder to save results
    mkdir(r_folder)

    datainfo=dir(path{p});
    fname{p}={datainfo.name};
    fname{p}(~contains(fname{p},'.tif'))=[];

    exps{p}=cellfun(@(x) x(1), fname{p});
    exp_labels=unique(exps{p});
    for e=1:numel(exp_labels)
        nexp{p}(e)=sum(exps{p}==exp_labels(e));
    end
    nfile{p}=numel(fname{p});

    for i=1:nfile{p}
       colocal{p}{i}=[]; 
       spot_I{p}{i}=[]; mean_I{p}{i}=[];

       img=strcat(path{p},'\',fname{p}{i});
       info=imfinfo(img);

       %image loading
       for imageNumber=1:size(info,1)
           img_raw(:,:,imageNumber)=imread(img,'index',imageNumber);
       end

       if i==find((exps{p}==fname{p}{i}(1)),1)
           background=bkg;
           leakage=0.15;
%            [tform,comp_img]=find_tform(img_raw(:,:,1)-background,disp);
       end

       %colocalization
       [up,down]=dualviewer_spliter(img_raw-background,tform);   
       down=down-up*leakage; %leakage substraction
       [co_pos{p}{i},~]=colocalization(strcat(r_folder,'\colocal_',fname{p}{i},'.png'),up(r0:r0+nr-1,c0:c0+nc-1,target_frame), down(r0:r0+nr-1,c0:c0+nc-1,target_frame),thres,max_off,disp);

       writeTIFrgb(down,up,zeros(size(up)),strcat(r_folder,'\composite_',fname{p}{i}));
        
       %% spot intensity analysis
       [nmol{p}{i},xpos,ypos,~,~]=particle_detection(up(r0:r0+nr-1,c0:c0+nc-1,target_frame),thres(1),0);
       mean_I{p}{i}=mean(up(r0:r0+nr-1,c0:c0+nc-1,target_frame),"all");
       for mol=1:nmol{p}{i}
           r=round(ypos(mol)+(-2:2)); c=round(xpos(mol)+(-2:2));              
           spot_I{p}{i}(mol,:)=[mean(up(r0+r-1,c0+c-1,target_frame),'all'),mean(down(r0+r-1,c0+c-1,target_frame),'all')];
           if isempty(cell2mat(co_pos{p}{i}{1})) colocal{p}{i}(mol)=0; 
           else colocal{p}{i}(mol)=sum(ismember([xpos(mol) ypos(mol)],cell2mat(co_pos{p}{i}{1}),'rows'));
           end
       end

      %% spot intensity analysis - DID
       [nmol{p}{i},xpos,ypos,~,~]=particle_detection(down(r0:r0+nr-1,c0:c0+nc-1,target_frame),thres(2),1);
       mean_I2{p}{i}=mean(down(r0:r0+nr-1,c0:c0+nc-1,target_frame),"all");
       for mol=1:nmol{p}{i}
           r=round(ypos(mol)+(-2:2)); c=round(xpos(mol)+(-2:2));              
           spot_I2{p}{i}(mol,:)=[mean(up(r0+r-1,c0+c-1,target_frame),'all'),mean(down(r0+r-1,c0+c-1,target_frame),'all')];
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
figure(Position=[500 300 300 400]);

for p=1:numel(path)
    temp1{p}=[]; temp2{p}=[];
    for i=1:nfile{p}
        temp1{p}=[temp1{p}; nnz(colocal{p}{i})/numel(colocal{p}{i})];
        temp2{p}=[temp2{p}; mean_I{p}{i}];
    end    
end

plist=[4,3,1];
for e=1:3
    plotSpread(temp1{plist(e)},'xValues',e, 'SpreadWidth',0.4, 'distributionColors','k')
    boxplot(temp1{plist(e)},'Positions',e,'Symbol','');
end
ylim([0,1])
ylabel('VAMP to ΔN ratio')

figure(Position=[500 300 300 400]);
for e=1:3
    plotSpread(temp2{plist(e)},'xValues',e, 'SpreadWidth',0.4, 'distributionColors','k')
    boxplot(temp2{plist(e)},'Positions',e,'Symbol','');
end
ylabel('mean intensity')
ylim([0,250])
