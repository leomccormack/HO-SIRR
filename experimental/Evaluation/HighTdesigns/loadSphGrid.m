function grid = loadSphGrid(filepath)
%LOADSPHGRID Summary of this function goes here
%   Detailed explanation goes here

data = importdata(filepath);

vecs = data(3:end);
grid.nPoints = length(vecs)/3;
vecs = reshape(vecs,3,grid.nPoints).';
grid.aziElev = unitCart2sph(vecs);
grid.xyz = vecs;

grid.faces = sphDelaunay(grid.aziElev);

% find average angular distance between adjacent points
avg_angle = 0;
nFaces = size(grid.faces,1);
for nf=1:nFaces
    face_idx = grid.faces(nf,:);
    angle1 = acos(dot(vecs(face_idx(1),:), vecs(face_idx(2),:)));
    angle2 = acos(dot(vecs(face_idx(2),:), vecs(face_idx(3),:)));
    angle3 = acos(dot(vecs(face_idx(3),:), vecs(face_idx(1),:)));
    avg_angle = avg_angle + (angle1+angle2+angle3)/3;
end
grid.avgAngle = avg_angle/nFaces;

end

