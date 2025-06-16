clear all; close all;
%% adjustable
r0=1; c0=200;
nr=200 ; nc=200 ;
sz=3;
%parameters for colocalization

target_frame=1;

thres=[15 20];
disp=[100, 700]; 

max_off=4; 

%%
plist=dir;
isSubfolder=[plist.isdir]; isSubfolder(1:2)=0;
path={plist(isSubfolder).name};
path(contains(path, 'result'))=[];
tforms{1}=load('tform1.mat').tform; tforms{2}=load('tform2.mat').tform;
% path={'raw2'};
for p=1:size(path,2)
    r_folder=['result_' path{p}]; %folder to save results
    mkdir(r_folder)

    datainfo=dir(path{p});
    fname{p}={datainfo.name};
    fname{p}(~contains(fname{p},'.tif'))=[];

    exps{p}=cellfun(@(x) x(1), fname{p});
    exp_labels{p}=unique(exps{p});
    for e=1:numel(exp_labels{p})
        nexp{p}(e)=sum(exps{p}==exp_labels{p}(e));
    end
    nfile{p}=numel(fname{p});

    %% background calculation 
    bkg={};
    for i=1:nfile{p}
       bkg_temp{p}{i}=[];
       img=strcat(path{p},'\',fname{p}{i});
       info=imfinfo(img);

       %image loading
       for imageNumber=1:size(info,1)
           img_raw(:,:,imageNumber)=imread(img,'index',imageNumber);
       end
       bkg_temp{p}{i}=min(img_raw,[],3);
    end

    for e=1:numel(exp_labels{p})
        bkg_stack=bkg_temp{p}(exps{p}==exp_labels{p}(e))';  
        bkg{p}{e}=min(cat(3,bkg_stack{:}),[],3); 
    end

    clear bkg_temp
    %%

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
           background=bkg{p}{find((exp_labels{p}==fname{p}{i}(1)),1)};
           leakage=0;
           tform=tforms{find((exp_labels{p}==fname{p}{i}(1)),1)};
%            [tform,comp_img]=find_tform(img_raw(:,:,1)-background,disp);
       end

       %colocalization
       [up,down]=dualviewer_spliter(img_raw-background,tform);   
       down=down-up*leakage; %leakage substraction
       [co_pos{p}{i},~]=colocalization(strcat(r_folder,'\colocal_',fname{p}{i},'.png'),up(r0:r0+nr-1,c0:c0+nc-1,target_frame), down(r0:r0+nr-1,c0:c0+nc-1,target_frame),thres,max_off,disp);

       writeTIFrgb(zeros(size(up)),up,down,strcat(r_folder,'\composite_',fname{p}{i}));
        
       %% spot intensity analysis
       [nmol{p}{i},xpos,ypos,~,~]=particle_detection(up(r0:r0+nr-1,c0:c0+nc-1,target_frame),thres(1),0);
       for mol=1:nmol{p}{i}
           r=round(ypos(mol)+(-sz:sz)); c=round(xpos(mol)+(-sz:sz));              
           spot_I{p}{i}(mol,:)=[mean(up(r0+r-1,c0+c-1,target_frame),'all'),mean(down(r0+r-1,c0+c-1,target_frame),'all')];
           if isempty(cell2mat(co_pos{p}{i}{1})) colocal{p}{i}(mol)=0; 
           else colocal{p}{i}(mol)=sum(ismember([xpos(mol) ypos(mol)],cell2mat(co_pos{p}{i}{1}),'rows'));
           end
       end

       %%
       clear img_raw up down mov
    end
end

save('data.mat')


%% data loading
clear
load data.mat
%%
figure;
edges=0:80:1600;

tmpdat{1}=cell2mat(spot_I{1}(1:5)'); tmpdat{2}=cell2mat(spot_I{1}(6:10)');
histogram(tmpdat{2}(:,1),edges,'Normalization','probability','EdgeAlpha',0); hold on;
histogram(tmpdat{1}(:,1),edges,'Normalization','probability','EdgeAlpha',0); hold off;
legend({'-Syn','+Syn'})
ylabel('probability')
xlabel('intensity (a.u.)')

%%
figure;
edges=0:60:1600;

tmpdat{1}=cell2mat(spot_I{1}(1:5)'); tmpdat{2}=cell2mat(spot_I{1}(6:10)');
subplot(2,1,1)
histogram(tmpdat{2}(:,1),edges,'EdgeAlpha',0); hold on;
subplot(2,1,2)
histogram(tmpdat{1}(:,1),edges,'EdgeAlpha',0); hold off;

ylabel('count')
xlabel('intensity (a.u.)')
%%
tmp=[];
for i=1:5
    tmp=vertcat(tmp,sum(colocal{1}{i})/numel(colocal{1}{i}));
end
tmpdat{1}=tmp;

tmp=[];
for i=6:10
    tmp=vertcat(tmp,sum(colocal{1}{i})/numel(colocal{1}{i}));
end
tmpdat{2}=tmp;

figure;
for e=[2 1]
    bar(3-e,mean(tmpdat{e})); hold on
    errorbar(3-e,mean(tmpdat{e}),std(tmpdat{e}),'k')
end
ylim([0 1]);