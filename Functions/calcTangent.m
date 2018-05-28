function [ coords_Tangent ] = calcTangent( hullcortex, c, coords, dims, lngth,hemi)
%This function calculates the tangent plane from a point c and rotates the
%input matrix (existing out of (X,Y,Z) coordinates which specify the grid in a fixed 2D plane)
% so that it lies within the tangent plane to the cortex specified. 
% Also provive the dimension dims and the length of mesh that you input.

N = normElec(hullcortex,c,25);
N = N(1,:);

%calculate plane equation
d=dot(N,c);

point = c;
normal = N;

yy = reshape(coords(:,2),dims(1,1),dims(1,2));
zz = reshape(coords(:,3),dims(1,1),dims(1,2));

x=(-normal(2)*yy-normal(3)*zz -d)/normal(1);


tangPlane = [];
tangPlane(:,1)=reshape(x,lngth,1);
tangPlane(:,2)=reshape(yy,lngth,1);
tangPlane(:,3)=reshape(zz,lngth,1);


%Calculate the distance vectors M and N between two points in the different planes
N = tangPlane(1,:) - tangPlane(8,:) ;
M = coords(1,:) - coords(8,:);

%Calculate the angle between the planes 
AngleBetweenPlanes = atan2d(norm(cross(M,N)),dot(M,N));

%use the function [Point,IntersectionPlanes,~] = plane_intersect(...) to
%calculate the direction vector and to calculate a point on the straight line
%intersection of the planes
[Point,IntersectionPlanes,~] = plane_intersect([1 0 0],coords(1,:),normal,point);

%Create rotationmatrix with AxelRot function
RotMatrix = AxelRot(AngleBetweenPlanes,IntersectionPlanes,Point);


%Rotate the electrode locations from one plane to the perpendicular plane
BluePlane = coords;
RedPlane = [];
for l=1:length(BluePlane)
    auxRedPlane = RotMatrix * [BluePlane(l,:)'; 1];
    RedPlane(l,:) = auxRedPlane(1:3);
end

coords_Tangent = RedPlane;

if hemi == 'r'
    coorX = reshape(coords_Tangent(:,1),dims(1,1),dims(1,2));
    coorY = reshape(coords_Tangent(:,2),dims(1,1),dims(1,2));
    coorZ = reshape(coords_Tangent(:,3),dims(1,1),dims(1,2));

    coorX2=rot90(coorX,2);

    tangPlaneAlt(:,1)=reshape(coorX2,lngth,1);
    tangPlaneAlt(:,2)=reshape(coorY,lngth,1);
    tangPlaneAlt(:,3)=reshape(coorZ,lngth,1);
    coords_Tangent = tangPlaneAlt;
end
end

