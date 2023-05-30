%%%%%%%%%%% Load Data %%%%%%%%%%%
clear all
datapath = 'roi\dHCP_pipeline\fetal_atlas_changed_code_zrk\fetal\data\derivatives\';
prefix = 'ses-session1\anat\Native\';
CRLpath = 'fetus_segmentation_CRL\data\subject\';
SubjList=dir(CRLpath);
N = length(SubjList);
resolution = 0.8;
for i = 3:N
    s=char(SubjList(i).name);
    week=str2double(s(1:2));
    day=str2double(s(4));
    week_day=week+day/7;
    SubjList(i).GA=week_day; 
end

%%%%%%%%%%% Regional %%%%%%%%%%%
load Region_CRL

for i = 3:N
        % load surface files
        WM_L = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_left_white.surf.gii']);
        Sulc_L = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_left_sulc.shape.gii']);
        MeanCur_L = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_left_mean_curvature.shape.gii']);
        GaussCur_L = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_left_gauss_curvature.shape.gii']);
        Curness_L = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_left_curvedness.shape.gii']);
        K1_L = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_left_K1.shape.gii']);
        K2_L = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_left_K2.shape.gii']);
        Thick_L = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_left_thickness.shape.gii']);
        WM_R = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_right_white.surf.gii']);
        Sulc_R = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_right_sulc.shape.gii']);
        MeanCur_R = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_right_mean_curvature.shape.gii']);
        GaussCur_R = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_right_gauss_curvature.shape.gii']);
        Curness_R = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_right_curvedness.shape.gii']);
        K1_R = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_right_K1.shape.gii']);
        K2_R = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_right_K2.shape.gii']);
        Thick_R = gifti([datapath,'sub-',SubjList(i).name,'\',prefix,'sub-',SubjList(i).name,'_ses-session1_right_thickness.shape.gii']);
        Regional_L = gifti([CRLpath,SubjList(i).name,'\surface\sub-',SubjList(i).name,'_ses-session1_left_ROI_dilD_merge.label.gii']);
        Regional_R = gifti([CRLpath,SubjList(i).name,'\surface\sub-',SubjList(i).name,'_ses-session1_right_ROI_dilD_merge.label.gii']);
        SurfArea_L = gifti([CRLpath,SubjList(i).name,'\surface\sub-',SubjList(i).name,'_ses-session1_left_surface_area.shape.gii']);
        SurfArea_R = gifti([CRLpath,SubjList(i).name,'\surface\sub-',SubjList(i).name,'_ses-session1_right_surface_area.shape.gii']);
        
        % load volume file
        AllLabels = load_untouch_nii([CRLpath,SubjList(i).name,'\','segmentation','\',SubjList(i).name,'_cGM.nii.gz']);
        
        for k = 1:2:78
                % surface measure
                Region(k).index = find(Regional_L.cdata==Region(k).labelkey);
                Region(k).NumberVertices = length(Region(k).index);
                Region(k).Thickness = sum(Thick_L.cdata(Region(k).index))/Region(k).NumberVertices;
                Region(k).SulcalDepth = sum(abs(Sulc_L.cdata(Region(k).index)))/Region(k).NumberVertices;
                Region(k).MeanCurvature = sum(abs(MeanCur_L.cdata(Region(k).index)))/Region(k).NumberVertices;
                Region(k).GaussCurvature = sum(abs(GaussCur_L.cdata(Region(k).index)))/Region(k).NumberVertices;
                Region(k).Curvedness = sum(abs(Curness_L.cdata(Region(k).index)))/Region(k).NumberVertices;
                Region(k).SurfaceArea = sum(SurfArea_L.cdata(Region(k).index));
        end
        
        for k = 2:2:78
                Region(k).index = find(Regional_R.cdata==Region(k).labelkey);
                Region(k).NumberVertices = length(Region(k).index);
                Region(k).Thickness = sum(Thick_R.cdata(Region(k).index))/Region(k).NumberVertices;
                Region(k).SulcalDepth = sum(abs(Sulc_R.cdata(Region(k).index)))/Region(k).NumberVertices;
                Region(k).MeanCurvature = sum(abs(MeanCur_R.cdata(Region(k).index)))/Region(k).NumberVertices;
                Region(k).GaussCurvature = sum(abs(GaussCur_R.cdata(Region(k).index)))/Region(k).NumberVertices;
                Region(k).Curvedness = sum(abs(Curness_R.cdata(Region(k).index)))/Region(k).NumberVertices;
                Region(k).SurfaceArea = sum(SurfArea_R.cdata(Region(k).index));
        end
        
        for k=1:78
            % volume
            logi = (AllLabels.img==Region(k).labelkey);
            Region(k).Volume = sum(logi(:))*(0.8^3);
        end
        
        % save data in SubjList
        for ii = 1:78
            SubjList(i).(['NV_',cell2mat(Region(ii).name)]) = Region(ii).NumberVertices;
        end
        
        for ii = 1:78
            SubjList(i).(['T_',cell2mat(Region(ii).name)]) = Region(ii).Thickness;
        end
        
        for ii = 1:78
            SubjList(i).(['S_',cell2mat(Region(ii).name)]) = Region(ii).SulcalDepth;
        end
        
        for ii = 1:78
            SubjList(i).(['M_',cell2mat(Region(ii).name)]) = Region(ii).MeanCurvature;
        end
        
        for ii = 1:78
            SubjList(i).(['G_',cell2mat(Region(ii).name)]) = Region(ii).GaussCurvature;
        end
        
        for ii = 1:78
            SubjList(i).(['C_',cell2mat(Region(ii).name)]) = Region(ii).Curvedness;
        end
        
        for ii = 1:78
            SubjList(i).(['V_',cell2mat(Region(ii).name)]) = Region(ii).Volume;
        end
        
        for ii = 1:78
            SubjList(i).(['SA_',cell2mat(Region(ii).name)]) = Region(ii).SurfaceArea;
        end
end

%% save result
cd('CorticalAnalysis\')
% delete result.xlsx
fieldname=fields(SubjList);
writecell(fieldname','result(1).xlsx')
writecell((struct2cell(SubjList(3:(N)))'),'result(1).xlsx','WriteMode','append')
