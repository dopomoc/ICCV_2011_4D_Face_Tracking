
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% paint_cylindrical_UVs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script takes as input:
% - a sequence of OBJ files with a stereo UV mapping
% - sequences of 1C/2C colour stereo camera images 
% - camera calibration data
% - cylindrical unwrapping data (new OBJ file/UV coordinates)
%
% This script outputs:
%
% - a new set of cylindrical UV maps created from the 1C/2C views and
% stereo OBJ file. These have file name AU_[NAME]_[IM_NUMBER]_w.bmp
%
% The overall purpose of this script is to take raw data from the cameras -
% OBJs and images, and create UV maps with a cylindrical UV mapping. The UV
% coords for the cylindrical mapping are known, so the idea is to select
% the best faces to select from 1C or 2C, and paint them onto the
% cylindrical UV map in the right place.
%
% One reason for this script is that the original stereo UVs are a little
% noisy, so the cylindrical unwrapping also creates noisy UVs. This script
% creates high quality UVs from the original 1C/2C camera images using the
% camera calibration matrix for the projection.
%
% Method:
%
% - Raster scan unwrapped UV coordinates.
% - For each pixel inside bounding box of u,v map:
%    -what triangle (unwrapped uv) does it belong to?
%    -what is the barycentric coordinate (unwrapped uv)?
%    -what does the pixel look like at this barycentric coordinate in the
%    stereo uv map (inside matching triangle in stereo map).
%    -write this pixel value to unwrapped uv map
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIX_JOIN=1; % Smooth the join in the new cylindrical UV map?
FIX_METHOD = 'AREA_AND_BOUND'; % Which method to select faces from 1C/2C?
FIX_WIDTH=5; % Smoothness filter width
newUVr = 1024; % Size of new UV map (rows)
newUVc = 1280; % Size of new UV map (columns)
FILTER_ITERS = 3; % Number of smoothing iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load camera calibration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NB - hard coded this for demo purposes - you could also load these in

% Load 1C calib matrix
Cal_1C.M = [0.858401484 0.228175161 -0.459437687 ; ...
    -0.0492414924 -0.854840292 -0.516549465; ...
    -0.510609604 0.466030225 -0.722560490 ]';
Cal_1C.X = 561.595017;
Cal_1C.Y =  -488.097124;
Cal_1C.Z = 894.552876;
Cal_1C.f = 20.9742442;
Cal_1C.K = 0.000444987580;
Cal_1C.K2 = -1.15206946e-005;
Cal_1C.S = 1.0;
Cal_1C.x_ = 0.00465000000;
Cal_1C.y_ = 0.00465000000;
Cal_1C.a = 661.532623;
Cal_1C.b = 484.318751;
Cal_1C.is = [1280 1024];
Cal_1C.c = 0.00000000;

% Load 2C calib matrix
Cal_2C.M = [0.774972732 -0.331350090 0.538167615  ; ...
    -0.0225409314 -0.865488103 -0.500422072 ; ...
    0.631592567 0.375682661 -0.678198620 ]';
Cal_2C.X = -709.212682;
Cal_2C.Y =  -398.935198;
Cal_2C.Z = 813.592684;
Cal_2C.f = 20.9110983;
Cal_2C.K = 0.000209355627;
Cal_2C.K2 = 3.78815746e-006;
Cal_2C.S = 1.0;
Cal_2C.x_ = 0.00465000000;
Cal_2C.y_ = 0.00465000000;
Cal_2C.a = 601.705107;
Cal_2C.b = 501.748982;
Cal_2C.is = [1280 1024];
Cal_2C.c = 0.00000000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load new unwrapped OBJ data (files names), and blank images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NB - hard coded some of this for demo purposes

pathname=uigetdir(pwd,'Select a directory');
contentsOBJ = dir(fullfile(pathname,'*_w.obj'));
contentsBMP = dir(fullfile(pathname,'*.bmp'));

contentsTXT_coords = dir(fullfile(pathname,'*w_tcoords.txt'));
contentsTXT_verts = dir(fullfile(pathname,'*w_verts.txt'));

% Load tracked points - this is a guide (used as a BOUND), roughly where
% the nose should be
load([pathname '\mid_OFLOW_tracked.mat']);

files_1C = dir([pathname '\texture_1C*.bmp']);
files_2C = dir([pathname '\texture_2C*.bmp']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop following algorithm above - 1 iteration per UV map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for main_loop=1:size(contentsBMP,1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load raw OBJ with stereo UV mapping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currentOBJFileName = ['\' contentsOBJ(main_loop).name];
    currentOBJFileName = [pathname currentOBJFileName]
    [pathstr,name,ext,versn] = fileparts(currentOBJFileName);

    fidClean=fopen([currentOBJFileName],'r');

    disp('Loading mesh faces');
    tic;
    uv_unwrapped=[];
    faces_unwrapped=[];
    counter=1;
    while 1
        try
            tline = fgetl(fidClean);
            if ~ischar(tline), break, end

            matches_3 = findstr(tline, 'f');
            num_3 = length(matches_3);

            if num_3 > 0 & counter>4
                nums = sscanf(tline,'%*s %d %[/] %d %[/] %d %d %[/] %d %[/] %d %d %[/] %d %[/] %d');
                nums = nums([1 6 11]);
                faces_unwrapped = [faces_unwrapped nums];
            end

            counter=counter+1;
        catch
            disp('bad symbol');
        end
    end
    fclose(fidClean);
    toc;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load stereo UV image and text coords (taken from OBJ)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    currentBMPFileName = ['\' contentsBMP(main_loop).name];
    currentBMPFileName = [pathname currentBMPFileName]
    [pathstr,name,ext,versn] = fileparts(currentBMPFileName);

    % Load stereo image
    stereo_image = imread([currentBMPFileName(1:end-4) '.bmp']);
    % Create new cylindrical unwrapped image of the same size
    unwrapped_image = zeros(size(stereo_image,1),size(stereo_image,2),size(stereo_image,3));
    
    [unwrapped_rows_orig,unwrapped_cols_orig,NN] = size(unwrapped_image);

    uv_unwrapped = load([pathname '\' contentsTXT_coords(main_loop).name])';
    uv_unwrapped(2,:) = (1-(uv_unwrapped(2,:)/size(unwrapped_image,1)))*size(unwrapped_image,1);

    unwrapped_image = imresize(unwrapped_image,[newUVr newUVc]);

    [unwrapped_rows,unwrapped_cols,NN] = size(unwrapped_image);
    [stereo_rows,stereo_cols,NN] = size(stereo_image);

    uv_unwrapped(1,:) = (uv_unwrapped(1,:)/unwrapped_cols_orig)*unwrapped_cols;
    uv_unwrapped(2,:) = (uv_unwrapped(2,:)/unwrapped_rows_orig)*unwrapped_rows;

    coords_unwrapped = uv_unwrapped;

    temp = contentsTXT_coords(main_loop).name;
    uv_stereo = load([pathname '\' temp(1:end-13) 'tcoords.txt'])';
    uv_stereo(2,:) = (1-(uv_stereo(2,:)/size(stereo_image,1)))*size(stereo_image,1);

    coords_stereo = uv_stereo;

    faces_stereo=faces_unwrapped;
 
    disp('Calculating 1C/2C UV coordinates..');
  
    % now get 1C and 2C coords for this
    Im_1C = imread([pathname '\' files_1C(main_loop).name]);
    Im_2C = imread([pathname '\' files_2C(main_loop).name]);

    [Im_1C_r,Im_1C_c,err] = size(Im_1C);

    temp = contentsTXT_verts(main_loop).name;
    clean_vertices = load([pathname '\' temp(1:end-12) '_verts.txt'])';

    % find 1C and 2D indices
    inds_1C = find(uv_stereo(1,:)<(0.49*size(stereo_image,2)));
    inds_2C = find(uv_stereo(1,:)>(0.48*size(stereo_image,2)));

    all_inds_1C{main_loop} = inds_1C;
    all_inds_2C{main_loop} = inds_2C;

    verts_1C = [clean_vertices(1,inds_1C);clean_vertices(2,inds_1C);clean_vertices(3,inds_1C)];
    verts_2C = [clean_vertices(1,inds_2C);clean_vertices(2,inds_2C);clean_vertices(3,inds_2C)];


    coords_unwrapped = [uv_unwrapped(1,:)*unwrapped_cols; (1 - uv_unwrapped(2,:))*unwrapped_rows];

    new_unwrapped = zeros(unwrapped_rows,unwrapped_cols,3);
    debug_stereo = zeros(stereo_rows,stereo_cols,1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each face in the new UV map, select 1C or 2C to paint from
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Writing face from 1C/2C UV to unwrapped UV..');
    tic

    Faces_in_1C = [];
    Faces_in_2C = [];

    for face_count = 1:1:size(faces_unwrapped,2)

        this_uv_unwrapped_face = [uv_unwrapped(:,faces_unwrapped(1,face_count));...
            uv_unwrapped(:,faces_unwrapped(2,face_count));...
            uv_unwrapped(:,faces_unwrapped(3,face_count))];

        this_face = [faces_unwrapped(1,face_count);faces_unwrapped(2,face_count);faces_unwrapped(3,face_count)];

        min_u = round(min(this_uv_unwrapped_face(1:2:end,1)));
        max_u = round(max(this_uv_unwrapped_face(1:2:end,1)));
        min_v = round(min(this_uv_unwrapped_face(2:2:end,1)));
        max_v = round(max(this_uv_unwrapped_face(2:2:end,1)));


        this_uv_stereo_face = [uv_stereo(:,faces_stereo(1,face_count));...
            uv_stereo(:,faces_stereo(2,face_count));...
            uv_stereo(:,faces_stereo(3,face_count))];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % There are lots of ways to paint the cylindrical UV map.
        % Basically, faces should be chosen from 1C or 2C and painted
        % on to the cylindrical UV map. How this is done is open, and one
        % method only is highlighted below as that's the method I finally
        % went with.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        switch FIX_METHOD
                
               
            case 'AREA_AND_BOUND'

                
                this_face_uv = [faces_unwrapped(1,face_count);faces_unwrapped(2,face_count);faces_unwrapped(3,face_count)];

                try
                    BOUND = (((rawTrainingData(1,main_loop)/720)+(rawTrainingData(3,main_loop)/720))/2)*unwrapped_cols;
                catch
                    BOUND = (rawTrainingData(1,main_loop)/720)*unwrapped_cols;
                end

                v1.x = clean_vertices(1,this_face_uv(1,1));
                v1.y = clean_vertices(2,this_face_uv(1,1));
                v1.z = clean_vertices(3,this_face_uv(1,1));

                v2.x = clean_vertices(1,this_face_uv(2,1));
                v2.y = clean_vertices(2,this_face_uv(2,1));
                v2.z = clean_vertices(3,this_face_uv(2,1));

                v3.x = clean_vertices(1,this_face_uv(3,1));
                v3.y = clean_vertices(2,this_face_uv(3,1));
                v3.z = clean_vertices(3,this_face_uv(3,1));

                % Use the tsai projection method to find where the face in
                % the cylindrical UV coords are in the 1C and 2C image. We
                % can then decide which image to take from - either 1C or
                % 2D - later on.
                p1_1C = tsai(Cal_1C,v1);
                p2_1C = tsai(Cal_1C,v2);
                p3_1C = tsai(Cal_1C,v3);

                p1_2C = tsai(Cal_2C,v1);
                p2_2C = tsai(Cal_2C,v2);
                p3_2C = tsai(Cal_2C,v3);

                this_uv_1C_face = [p1_1C.x;p1_1C.y;...
                    p2_1C.x;p2_1C.y;...
                    p3_1C.x;p3_1C.y];

                this_uv_2C_face = [p1_2C.x;p1_2C.y;...
                    p2_2C.x;p2_2C.y;...
                    p3_2C.x;p3_2C.y];

                this_uv_1C_face(1:2:end,:) = (this_uv_1C_face(1:2:end,:)/1280)*Im_1C_c;
                this_uv_1C_face(2:2:end,:) = (1-(this_uv_1C_face(2:2:end,:)/1024))*Im_1C_r;

                this_uv_2C_face(1:2:end,:) = (this_uv_2C_face(1:2:end,:)/1280)*Im_1C_c;
                this_uv_2C_face(2:2:end,:) = (1-(this_uv_2C_face(2:2:end,:)/1024))*Im_1C_r;

                % Get the area of the face in 1C and 2C
                Area_1C = triangle_area([p1_1C.x p1_1C.y p1_1C.z; p2_1C.x p2_1C.y p2_1C.z; p3_1C.x p3_1C.y p3_1C.z]);
                Area_2C = triangle_area([p1_2C.x p1_2C.y p1_2C.z; p2_2C.x p2_2C.y p2_2C.z; p3_2C.x p3_2C.y p3_2C.z]);

                % Does the face look toward the camera? Cross produce face
                % normal with camera direction to check.
                cp_1C = cross([p1_1C.x;p1_1C.y;p1_1C.z],[p2_1C.x;p2_1C.y;p2_1C.z]);
                cp_2C = cross([p1_2C.x;p1_2C.y;p1_2C.z],[p2_2C.x;p2_2C.y;p2_2C.z]);
                dp_1C = min(dot(cp_1C,[Cal_1C.X;Cal_1C.Y;Cal_1C.Z]),dot(cp_1C,-[Cal_1C.X;Cal_1C.Y;Cal_1C.Z]));
                dp_2C = min(dot(cp_2C,[Cal_2C.X;Cal_2C.Y;Cal_2C.Z]),dot(cp_2C,-[Cal_2C.X;Cal_2C.Y;Cal_2C.Z]));


                this_uv_stereo_face = [uv_stereo(:,faces_stereo(1,face_count));...
                    uv_stereo(:,faces_stereo(2,face_count));...
                    uv_stereo(:,faces_stereo(3,face_count))];

              
                % Do checks - based on area and whether the face is above
                % or below the BOUND. Decide whether to take a face from 1C
                % or 2C, and use this to paint the unwrapped UV map.

                if(Area_1C>Area_2C) & mean(this_uv_unwrapped_face(1:2:end,1))<BOUND 
                                        
                    % It's in 2C
                    Faces_in_2C = [Faces_in_2C this_uv_unwrapped_face];
                    
                    try
                        for u = min_u:max_u
                            for v = min_v:max_v
                                in = inpolygon(u,v,this_uv_unwrapped_face(1:2:end,1),this_uv_unwrapped_face(2:2:end,1));

                                [ii,vv] = find(in>0);

                                if(ii>0)

                                    % what are barycentric coords for this coord?
                                    new_point = mapPointBetweenTriangles2D([u;v],reshape(this_uv_unwrapped_face,2,3)',reshape(this_uv_2C_face,2,3)');

                                    new_unwrapped(v,u,:) = Im_2C(round(new_point(2,1)),round(new_point(1,1)),:);

                                end
                            end
                        end
                    catch

                    end
                                        
                    
                elseif(Area_1C<Area_2C) & mean(this_uv_unwrapped_face(1:2:end,1))>BOUND

                    % It's in 1C
                    Faces_in_1C = [Faces_in_1C this_uv_unwrapped_face];
                    
                    try
                        for u = min_u:max_u
                            for v = min_v:max_v
                                in = inpolygon(u,v,this_uv_unwrapped_face(1:2:end,1),this_uv_unwrapped_face(2:2:end,1));

                                [ii,vv] = find(in>0);

                                if(ii>0)

                                    % what are barycentric coords for this coord?
                                    new_point = mapPointBetweenTriangles2D([u;v],reshape(this_uv_unwrapped_face,2,3)',reshape(this_uv_1C_face,2,3)');

                                    new_unwrapped(v,u,:) = Im_1C(round(new_point(2,1)),round(new_point(1,1)),:);

                                end
                            end
                        end
                    catch

                    end
                    
                    
                elseif(Area_1C>Area_2C) & mean(this_uv_unwrapped_face(1:2:end,1))>BOUND
                    
                    % It's in 1C
                    Faces_in_1C = [Faces_in_1C this_uv_unwrapped_face];
                    
                    try
                        for u = min_u:max_u
                            for v = min_v:max_v
                                in = inpolygon(u,v,this_uv_unwrapped_face(1:2:end,1),this_uv_unwrapped_face(2:2:end,1));

                                [ii,vv] = find(in>0);

                                if(ii>0)

                                    % what are barycentric coords for this coord?
                                    new_point = mapPointBetweenTriangles2D([u;v],reshape(this_uv_unwrapped_face,2,3)',reshape(this_uv_1C_face,2,3)');

                                    new_unwrapped(v,u,:) = Im_1C(round(new_point(2,1)),round(new_point(1,1)),:);

                                end
                            end
                        end
                    catch

                    end

                elseif(Area_1C<Area_2C) & mean(this_uv_unwrapped_face(1:2:end,1))<BOUND
                    
                    % It's in 2C
                    Faces_in_2C = [Faces_in_2C this_uv_unwrapped_face];
                    
                    try
                        for u = min_u:max_u
                            for v = min_v:max_v
                                in = inpolygon(u,v,this_uv_unwrapped_face(1:2:end,1),this_uv_unwrapped_face(2:2:end,1));

                                [ii,vv] = find(in>0);

                                if(ii>0)

                                    % what are barycentric coords for this coord?
                                    new_point = mapPointBetweenTriangles2D([u;v],reshape(this_uv_unwrapped_face,2,3)',reshape(this_uv_2C_face,2,3)');

                                    new_unwrapped(v,u,:) = Im_2C(round(new_point(2,1)),round(new_point(1,1)),:);

                                end
                            end
                        end
                    catch

                    end

                else

                    % Face is undecided - where to place..?
                    % This is more for debugging, as the cases above are
                    % always satisfied. 
                end
               
        end        

    end
    toc;

    
    
    disp('Smoothing UV join artefact..');
    tic;

    % Fix Join. Find the middle join (1C/2C faces that share vertices in the new cylindrical UV map).
    % Then, filter down this join.
    if(FIX_JOIN)

        Faces_in_1C_Mod = round(Faces_in_1C);
        Faces_in_2C_Mod = round(Faces_in_2C);

        temp_1C = reshape(Faces_in_1C_Mod,2,size(Faces_in_1C_Mod,2)*3);
        temp_2C = reshape(Faces_in_2C_Mod,2,size(Faces_in_2C_Mod,2)*3);

        matches=[];

        for i = 1:size(temp_1C,2)

            p1 = temp_1C(1:2,i);

            % is current 1C point in Faces_2C
            [i1,v1] = find(temp_2C(1,:)==p1(1,1));
            [i2,v2] = find(temp_2C(2,:)==p1(2,1));

            a = intersect(v1,v2);

            if(~isempty(a))
                matches=[matches temp_2C(:,a)];
            end

        end

        [ii,vv] = sort(matches(2,:));

        sorted_matches = matches(:,vv);


        % alter negative values if they exist
        for i=1:size(sorted_matches,2)
            if(sorted_matches(1,i)<1)
                sorted_matches(1,i)=1;
            end
            if(sorted_matches(2,i)<1)
                sorted_matches(2,i)=1;
            end
        end

        min_y = min(sorted_matches(2,:));
        max_y = max(sorted_matches(2,:));

        all_ys = [min_y:max_y];

        all_join_coords=zeros(2,1024);

        for i=1:size(sorted_matches,2)-1

            start_y = sorted_matches(2,i);
            end_y = sorted_matches(2,i+1);

            len = end_y-start_y;

            if(len>1)
         
                diff = (sorted_matches(1,i+1)-sorted_matches(1,i))/len;

                curr_x = sorted_matches(1,i);

                for j = start_y:end_y+1

                    all_join_coords(:,j) = [round(curr_x);j];
                    curr_x=curr_x+diff;

                end

            end

        end

        % Go down join and filter

        temp = new_unwrapped;
        

        for i=1:size(all_join_coords,2)

            if(all_join_coords(1,i)~=0)
                
                for filter_loop=1:FILTER_ITERS
                    
                    this_strip_R = temp(all_join_coords(2,i),all_join_coords(1,i)-FIX_WIDTH:all_join_coords(1,i)+FIX_WIDTH,1);
                    this_strip_G = temp(all_join_coords(2,i),all_join_coords(1,i)-FIX_WIDTH:all_join_coords(1,i)+FIX_WIDTH,2);
                    this_strip_B = temp(all_join_coords(2,i),all_join_coords(1,i)-FIX_WIDTH:all_join_coords(1,i)+FIX_WIDTH,3);
                    
                    temp(all_join_coords(2,i),all_join_coords(1,i)-FIX_WIDTH:all_join_coords(1,i)+FIX_WIDTH,1) = averageFilter1D(this_strip_R')';
                    temp(all_join_coords(2,i),all_join_coords(1,i)-FIX_WIDTH:all_join_coords(1,i)+FIX_WIDTH,2) = averageFilter1D(this_strip_G')';
                    temp(all_join_coords(2,i),all_join_coords(1,i)-FIX_WIDTH:all_join_coords(1,i)+FIX_WIDTH,3) = averageFilter1D(this_strip_B')';
                                        
                end
                
            end

        end

        toc;

        new_unwrapped = imresize(temp,[unwrapped_rows_orig unwrapped_cols_orig]);        
        imwrite(new_unwrapped/255,[currentBMPFileName(1:end-4) '_w.bmp'],'BMP');

    else

        new_unwrapped = imresize(new_unwrapped,[unwrapped_rows_orig unwrapped_cols_orig]);        
        imwrite(new_unwrapped/255,[currentBMPFileName(1:end-4) '_w.bmp'],'BMP');

    end

end




