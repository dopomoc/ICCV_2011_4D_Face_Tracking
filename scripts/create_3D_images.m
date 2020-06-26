

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create_3D_images 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script takes as input:
% - a sequence of OBJ files _w (cylindrical mapping)
% - a sequence of BMPs _w (cylindrical UV maps)
%
% This script outputs:
%
% - 3D images, i.e. images where each x,y is a x,y,z value calculated from
% the OBJ.
%
% - This 3D image format is pretty convenient. What it allows is for
% non-rigid alignment to be calculated based on the RGB values of the UV
% maps - say using optical flow, or AAM tracking (see ICCV paper), or other tracking - 
% and then the alignment data to be applied to the 3D image. Applying in
% this sense means alignment to a reference image - say a neutral
% expression. Once this is done, a canonical vertex mapping - or a single
% stable known facial vertex configuration - can be applied to the whole
% sequence of OBJs. Originally, the OBJS have different vertex topologies
% etc, so this fixes the problem. The quality of the topology and it's
% stability is directly proportional to the quality of the non-rigid
% alignment as applied to the RGB UV maps. So a better method can be
% plugged in at a later date. This script however, just creates the 3D
% images, which are used in non-rigid alignment in the next part of the
% demo.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathname2=uigetdir(pwd,'Select a directory of OBJs and UVs (cylindrical mapping).');

allAlphaVals = [];
allBetaVals = [];

contentsTXT_coords = dir(fullfile(pathname2,'*w_tcoords.txt'));
contentsTXT_verts = dir(fullfile(pathname2,'*w_verts.txt'));

contentsBMP = dir(fullfile(pathname2,'*_w.bmp'));
contentsOBJ = dir(fullfile(pathname2,'*_w.obj'));
   
% Set the size of the 3D Images. For speed, this can be made smaller (done
% here!)
dim_r = 512;%1024;
dim_c = 640;%1280;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop - for each OBJ/UV, create a 3D Image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(contentsOBJ,1)
           
    tic;
    
    % Load UV map, XYZ coords, and UV coords
    I = imread([pathname2 contentsBMP(i).name]);
    
    I = imresize(I,[dim_r dim_c]);
    
    I_orig = imread([pathname2 contentsBMP(i).name]);
  
    I_XYZ = load([pathname2 contentsTXT_verts(i).name])';
    
    I_UV = load([pathname2 contentsTXT_coords(i).name])';
    
    I_UV(2,:) = (1-(I_UV(2,:)/size(I_orig,1)))*size(I_orig,1);
    I_UV(1,:) = (I_UV(1,:)/size(I_orig,2))*dim_c;
    I_UV(2,:) = (I_UV(2,:)/size(I_orig,1))*dim_r;
        
    % Load faces from OBJs
    disp('Loading faces..');
    I_Faces=[];
    fidClean=fopen([pathname2 contentsOBJ(i).name],'r');
    counter=1;
    while 1
        tline = fgetl(fidClean);
        if ~ischar(tline), break, end
        
        matches_3 = findstr(tline, 'f');
        num_3 = length(matches_3);
        
        if num_3 > 0 & counter>4
            nums = sscanf(tline,'%*s %d %[/] %d %[/] %d %d %[/] %d %[/] %d %d %[/] %d %[/] %d');
            nums = nums([1 6 11]);
            I_Faces = [I_Faces nums];
        end
        counter=counter+1;
    end
    fclose(fidClean);
    
    I_3D = zeros(size(I,1),size(I,2),3);
    
    Coords_3D = [];
    
    % Go through the OBJ face list, use barycentric coordinates to estimate 3D XYZ
    % points using UV coords for a face, and 3D coordinates associated with
    % each face.
    disp('Going through each face and creating 3D points..');
    for ii_ = 1:size(I_Faces,2)
        
        this_face = I_Faces(:,ii_);
        
        this_face_uv = [I_UV(1:2,this_face(1,1));...
            I_UV(1:2,this_face(2,1));...
            I_UV(1:2,this_face(3,1))];        
        
        this_face_vertices = [I_XYZ(:,I_Faces(1,ii_));...
            I_XYZ(:,I_Faces(2,ii_)); ...
            I_XYZ(:,I_Faces(3,ii_))];
        
        % make contrived points in and around current face
        
        min_u = round(min(this_face_uv(1:2:end,1)));
        max_u = round(max(this_face_uv(1:2:end,1)));
        min_v = round(min(this_face_uv(2:2:end,1)));
        max_v = round(max(this_face_uv(2:2:end,1)));
        
        
        
        for u = min_u:2:max_u
            for v = min_v:2:max_v
                
                try
                    in = inpolygon(u,v,this_face_uv(1:2:end,1),this_face_uv(2:2:end,1));
                    
                    [ii,vv] = find(in>0);
                    
                    if(ii>0)
                        new_point = mapPointBetweenTriangles3D([u;v], ...
                            [this_face_uv(1:2:end,1) this_face_uv(2:2:end,1)],...
                            [this_face_vertices(1:3:end,1) this_face_vertices(2:3:end,1) this_face_vertices(3:3:end,1)]);
                        
                        I_3D(v,u,:) = new_point;
                        
                        Coords_3D = [Coords_3D new_point];
                        
                    end
                catch
                    % Catch a mapping outside the image - this is fine, so
                    % continue
                end
            end
        end
        
    end
    
    toc;
    
    disp('Interpolating missing points..');
    
    % Interpolate any blank coordinates remaining in I_3D
    
    min_r = ceil(min(I_UV(2,:)))+1;
    max_r = floor(max(I_UV(2,:)))-1;
    
    min_c = ceil(min(I_UV(1,:)))+1;
    max_c = floor(max(I_UV(1,:)))-1;
        
    for r = min_r:max_r/2        
        for c = min_c:max_c/2
            
            try                
                if(I_3D(r,c,1)==0)
                    
                    rm1 = I_3D(r-1,c,1);
                    rp1 = I_3D(r+1,c,1);
                    cm1 = I_3D(r,c-1,1);
                    cp1 = I_3D(r,c+1,1);
                    
                    all_vals = [rm1;rp1;cm1;cp1];
                    [ii,vv] = find(all_vals==0);
                    
                    if(size(ii,1)<4)
                        I_3D(r,c,1) = sum(all_vals)/(4-size(ii,1));
                        I_3D(r,c,2) = sum([I_3D(r-1,c,2);I_3D(r+1,c,2);I_3D(r,c+1,2);I_3D(r,c-1,2)])/(4-size(ii,1));
                        I_3D(r,c,3) = sum([I_3D(r-1,c,3);I_3D(r+1,c,3);I_3D(r,c+1,3);I_3D(r,c-1,3)])/(4-size(ii,1));
                    end
                end                
            catch
                % Catch a mapping outside the image - this is fine, so
                % continue
            end
        end
    end
    
    for r = max_r:-1:max_r/2
        for c = max_c:-1:max_c/2
            try
                if(I_3D(r,c,1)==0)
                    
                    rm1 = I_3D(r-1,c,1);
                    rp1 = I_3D(r+1,c,1);
                    cm1 = I_3D(r,c-1,1);
                    cp1 = I_3D(r,c+1,1);
                    
                    all_vals = [rm1;rp1;cm1;cp1];
                    [ii,vv] = find(all_vals==0);
                    
                    if(size(ii,1)<4)
                        I_3D(r,c,1) = sum(all_vals)/(4-size(ii,1));
                        I_3D(r,c,2) = sum([I_3D(r-1,c,2);I_3D(r+1,c,2);I_3D(r,c+1,2);I_3D(r,c-1,2)])/(4-size(ii,1));
                        I_3D(r,c,3) = sum([I_3D(r-1,c,3);I_3D(r+1,c,3);I_3D(r,c+1,3);I_3D(r,c-1,3)])/(4-size(ii,1));
                    end
                end
            catch
                % Catch a mapping outside the image - this is fine, so
                % continue
            end
        end
    end
    
    for r = max_r:-1:max_r/2
        for c = min_c:max_c/2
            try
                if(I_3D(r,c,1)==0)
                    
                    rm1 = I_3D(r-1,c,1);
                    rp1 = I_3D(r+1,c,1);
                    cm1 = I_3D(r,c-1,1);
                    cp1 = I_3D(r,c+1,1);
                    
                    all_vals = [rm1;rp1;cm1;cp1];
                    [ii,vv] = find(all_vals==0);
                    
                    if(size(ii,1)<4)
                        I_3D(r,c,1) = sum(all_vals)/(4-size(ii,1));
                        I_3D(r,c,2) = sum([I_3D(r-1,c,2);I_3D(r+1,c,2);I_3D(r,c+1,2);I_3D(r,c-1,2)])/(4-size(ii,1));
                        I_3D(r,c,3) = sum([I_3D(r-1,c,3);I_3D(r+1,c,3);I_3D(r,c+1,3);I_3D(r,c-1,3)])/(4-size(ii,1));
                    end
                end
            catch
                % Catch a mapping outside the image - this is fine, so
                % continue
            end
        end
    end
    
    for r = min_r:max_r/2
        for c = max_c:-1:max_c/2
            try
                if(I_3D(r,c,1)==0)
                    
                    rm1 = I_3D(r-1,c,1);
                    rp1 = I_3D(r+1,c,1);
                    cm1 = I_3D(r,c-1,1);
                    cp1 = I_3D(r,c+1,1);
                    
                    all_vals = [rm1;rp1;cm1;cp1];
                    [ii,vv] = find(all_vals==0);
                    
                    if(size(ii,1)<4)
                        I_3D(r,c,1) = sum(all_vals)/(4-size(ii,1));
                        I_3D(r,c,2) = sum([I_3D(r-1,c,2);I_3D(r+1,c,2);I_3D(r,c+1,2);I_3D(r,c-1,2)])/(4-size(ii,1));
                        I_3D(r,c,3) = sum([I_3D(r-1,c,3);I_3D(r+1,c,3);I_3D(r,c+1,3);I_3D(r,c-1,3)])/(4-size(ii,1));
                    end
                end
            catch
                % Catch a mapping outside the image - this is fine, so
                % continue
            end
        end
    end
    
    I_3D = real(I_3D);   
    temp_filename = [pathname2 contentsOBJ(i).name];
   
    save([temp_filename(1:end-4) '_3D_Image.mat'],'I_3D');        
    
end



