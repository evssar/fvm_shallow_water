clear all;
close all;
clc;


R=zeros(3,3);
L=zeros(3,3);


Lx=1;
Ly=1;
Nx=60;
Ny=60;
dx=(2*Lx)/Nx;
dy=(2*Ly)/Ny;
x = -1-dx:dx:1+dx; 
y = -1-dy:dy:1+dy; 
[X,Y]=meshgrid(x,y);

cells={};
k=1;
for i=1:(Nx-1)
    for j=1:(Ny-1)
        cell={};
        cell.r1(1)=X(i,j);
        cell.r1(2)=Y(i,j);
        cell.r1(3)=0;
        cell.r2(1)=X(i,j+1);
        cell.r2(2)=Y(i,j+1);
        cell.r2(3)=0;
        cell.r3(1)=X(i+1,j+1);
        cell.r3(2)=Y(i+1,j+1);
        cell.r3(3)=0;
        cell.r4(1)=X(i+1,j);
        cell.r4(2)=Y(i+1,j);
        cell.r4(3)=0;
        cell.r=1/4*(cell.r1+cell.r2+cell.r3+cell.r4);
        cell.d1=cell.r2-cell.r1;
        cell.d2=cell.r3-cell.r2;
        cell.d3=cell.r4-cell.r3;
        cell.d4=cell.r1-cell.r4;
        cell.max_d=max([norm(cell.d1),norm(cell.d2),norm(cell.d3),norm(cell.d4)]);
        cell.n=cross(cell.r3-cell.r1,cell.r4-cell.r2);
        cell.A=norm(cell.n);
        cell.n=cell.n/cell.A;
        cell.neighbours=zeros(1,4);
        cell.to_delete=0;
        cells{k}=cell;
        k=k+1;
        
    end
end

N_cells=length(cells);

for k=1:N_cells
   
    if cells{k}.r(1) >=-0.3 && cells{k}.r(1) <=-0.25 && abs(cells{k}.r(2)) > 0.2
        cells{k}.to_delete=1;
    end
    
end

N_faces=0;
N_vertices=0;

j=1;

for i=1:N_cells
    
    if cells{i}.to_delete==0
        
        vertices(4*j-3,:)=cells{i}.r1;
        vertices(4*j-2,:)=cells{i}.r2;
        vertices(4*j-1,:)=cells{i}.r3;
        vertices(4*j-0,:)=cells{i}.r4;
        
        
        faces(j,1)=4*j-3;
        faces(j,2)=4*j-2;
        faces(j,3)=4*j-1;
        faces(j,4)=4*j-0;
        
        types(j,1)=cells{i}.to_delete;
        
        j=j+1;
        
        
    end
end


figure(1);
patch('Vertices',vertices,'Faces',faces,'FaceColor','flat','FaceVertexCData',types);
grid on;
colorbar;
axis equal;


grid_file=fopen('dam.grid','w');

fprintf(grid_file,'%i,%i\n',length(vertices),length(faces));
for i=1:length(vertices)
   fprintf(grid_file,'%0.6f,%0.6f,%0.6f\n',vertices(i,1),vertices(i,2),vertices(i,3));
end
for i=1:length(faces)
    fprintf(grid_file,'%i,%i,%i,%i\n',4*i-4,4*i-3,4*i-2,4*i-1);
end
fclose(grid_file);

