

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit_canonical_mesh_to_raw_OBJs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script takes as input:
% - a sequence of 3D Images
% - a sequence of BMPs _w (cylindrical UV maps)
%
% This script outputs:
%
% - aligned UV maps (BMPS)
% - OBJ files with a single mesh tracked through the sequence
%
% - There are lots of ways to do the mesh fitting. Here, each UV map is
% warped to a reference (frame 1) using tracked points (AAM) and a Thin
% Plate Spline warp. You can also use optical flow, or any other tracker
% you like, but here I am using TPS as it's stable and simple.
% The warps are also applied to the 3D Images, and once they have been
% aligned they are then sampled using the topology required for the overall
% mesh. This can be artist based, but here Im just using a regular grid for
% simplicity - so one vertex for every other pixel. This can be made one
% vertex per pixel, or vertices wherever the artist likes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Sequence Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pathname2=uigetdir(pwd,'Select a directory of 3D Images and UV maps (cylindrical mapping).');
contentsBMP = dir(fullfile(pathname2,'*_w.bmp')); 
contents3D = dir(fullfile(pathname2,'*_w_3D_Image.mat'));
load([pathname2 '\tracked_face_points_576_720_res.mat']);

% The points were tracked on smaller res versions of the images, so Im just
% putting these limits in here.
dim_r = 576;
dim_c = 720;

% Warping parameters
interp.method = 'nearest'; % interpolation method
interp.radius = 5; % radius or median filter dimension
interp.power = 2; %power for inverse weighting interpolation method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop. Load UV map and 3D Image. Warp UV map to reference and apply
% same warp to the 3D Image.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for main_loop = 1:size(contentsBMP,1)
    
    
    % Load RGB image
    currentImageName = contentsBMP(main_loop).name
    currentImageName = [pathname2 '\' currentImageName];
    currentImage = imread(currentImageName);
    currentImage = imresize(currentImage,[dim_r dim_c]);
    
    currentShape = rawTrainingData(:,main_loop);
    
    orig_image = imread(currentImageName);
    
    rSignal = currentImage(:,:,1);
    gSignal = currentImage(:,:,2);
    bSignal = currentImage(:,:,3);
    
    % Load 3D image
    currentImageName = contents3D(main_loop).name
    currentImageName = [pathname2 '\' currentImageName];
    currentImage_3D = load(currentImageName);
    currentImage_3D = currentImage_3D.temp_3D;
    currentImage_3D = imresize(currentImage_3D,[dim_r dim_c]);
    
    x3D = currentImage_3D(:,:,1);
    y3D = currentImage_3D(:,:,2);
    z3D = currentImage_3D(:,:,3);
    
    border = [1 1 1 dim_r dim_c dim_r dim_c 1]';
   
    start = [currentShape;border];
    end_ = [reference_face_shape;border];
    
    %%{
    [X, Y, Xw, Yw, outH, outW, interp] = tpswarp_pre(rSignal',[size(rSignal,1) size(rSignal,2)],[start(1:2:end,1) start(2:2:end,1)],...
        [end_(1:2:end,1) end_(2:2:end,1)],interp);
    
    [rSignal,imgwr,map] = interp2d(X(:), Y(:), rSignal', Xw, Yw, outH, outW, interp);
    rSignal=rSignal';
    
    [gSignal,imgwr,map] = interp2d(X(:), Y(:), gSignal', Xw, Yw, outH, outW, interp);
    gSignal=gSignal';
    
    [bSignal,imgwr,map] = interp2d(X(:), Y(:), bSignal', Xw, Yw, outH, outW, interp);
    bSignal=bSignal';
    
    [x3D,imgwr,map] = interp2d(X(:), Y(:), x3D', Xw, Yw, outH, outW, interp);
    x3D=x3D';
    
    [y3D,imgwr,map] = interp2d(X(:), Y(:), y3D', Xw, Yw, outH, outW, interp);
    y3D=y3D';
    
    [z3D,imgwr,map] = interp2d(X(:), Y(:), z3D', Xw, Yw, outH, outW, interp);
    z3D=z3D';
    
    
    rgbImage = zeros(576,720,3);
    rgbImage(:,:,1) = rSignal;
    rgbImage(:,:,2) = gSignal;
    rgbImage(:,:,3) = bSignal;
    
    I_3D = zeros(576,720,3);
    I_3D(:,:,1) = x3D;
    I_3D(:,:,2) = y3D;
    I_3D(:,:,3) = z3D;
    
    % Prepare image for display/further processing
    rgbImage = double(rgbImage)/255;
    
    temp = [pathname2 '\' contentsBMP(main_loop).name];
    imwrite(rgbImage,[temp(1:end-4) '_aligned.bmp'],'bmp');
    
    temp = [pathname2 '\' contents3D(main_loop).name];
    save([temp(1:end-4) '_aligned.mat'],'I_3D');
    
end


% Get the names of the aligned 3D Images
files_3D_out = dir(fullfile(pathname2,['*_3D_Image_aligned.mat']));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create mesh structure.
% Here, we just use one vertex for every other pixel, in a regular grid.
% But the sampling can be changed to whatever you like..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STEP=2; % The step size. 1 = one vertex per pixel. It's just a regular sampling here.

min_u = 2;
max_u = size(rgbImage,2);
min_v = 2;
max_v = size(rgbImage,1);

temp_UV=zeros(2,size([min_u:STEP:max_u],2)+size([min_v:STEP:max_v],2));

count=1;

for ii=min_u:STEP:max_u
    for jj=min_v:STEP:max_v
        temp_UV(1,count) = ii;
        temp_UV(2,count) = jj;
        count=count+1;
    end
end

x_mean = (reference_face_shape(1:2:end,1)/720)*dim_c;
y_mean = (reference_face_shape(2:2:end,1)/567)*dim_r;

k = convhull(x_mean,y_mean);

in = inpolygon(temp_UV(1,:),temp_UV(2,:),x_mean(k),y_mean(k));

[ii,vv] = find(in==1);

x_in = temp_UV(1,vv);
y_in = temp_UV(2,vv);

I_new_UV = [x_in;y_in];

TRI = delaunay(I_new_UV(1,:),I_new_UV(2,:));
temp = TRI(:,1);
TRI(:,1) = TRI(:,3);
TRI(:,3) = temp;
I_new_Faces = TRI';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in the 3D data from the aligned 3D Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertices_3D = {};

for i = 1:size(files_3D_out,1)
    i
    I = load([pathname2 '\' files_3D_out(i).name]);
    I = I.I_3D;
    
    coords_3D = [];
    
    for j=1:size(I_new_UV,2)
        coord = [I(I_new_UV(2,j),I_new_UV(1,j),1);I(I_new_UV(2,j),I_new_UV(1,j),2);I(I_new_UV(2,j),I_new_UV(1,j),3)];
        coords_3D = [coords_3D coord];
    end
    
    vertices_3D{i} = coords_3D;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out new OBJ files, where each now has the same structure and the
% mesh should be tracked to the sequence. The limit of this accuracy is the
% tracking method - as previously explained. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(files_3D_out,1)
    
    
    % write new obj to disk using 'neutral' face and UV coord list
    temp = [pathname2 '\' contentsBMP(i).name];
    fidOBJ=fopen([temp(1:end-4) '_aligned.obj'],'w');
    
    
    for k=1:size(vertices_3D{i},2)
        fprintf(fidOBJ,'v %f %f %f\n',vertices_3D{i}(1,k),vertices_3D{i}(2,k),vertices_3D{i}(3,k));
    end
    
    for k=1:size(I_new_UV,2)
        fprintf(fidOBJ,'vt %f %f\n',I_new_UV(1,k)/720,(1-I_new_UV(2,k)/576));
    end
    
    for k=1:size(vertices_3D{i},2)
        fprintf(fidOBJ,'vn %f %f %f\n',vertices_3D{i}(1,k),vertices_3D{i}(2,k),vertices_3D{i}(3,k));
        %fprintf(fidOBJ,'vn %f %f %f\n',1,1,1);
    end
    
    for k=1:size(I_new_Faces,2)
        fprintf(fidOBJ,'f %i/%i/%i %i/%i/%i %i/%i/%i\n',...
            I_new_Faces(1,k),I_new_Faces(1,k),I_new_Faces(1,k),...
            I_new_Faces(2,k),I_new_Faces(2,k),I_new_Faces(2,k),...
            I_new_Faces(3,k),I_new_Faces(3,k),I_new_Faces(3,k));
    end

    
    fclose(fidOBJ);
    
end































