accuracy = [1, 5,0.922;
          1, 10, 0.916;
          1, 15, 0.927;
          1, 20, 0.932;
          1, 25, 0.933;
          1, 30, 0.938;
          1, 35, 0.934;
          1, 40, 0.944;
          1, 45, 0.944;
          1, 50, 0.946;

          2, 5, 0.924;
          2, 10,0.941;
          2, 15,0.942;
          2, 20,0.946;
          2, 25,0.938;
          2, 30,0.941;
          2, 35,0.943;
          2, 40,0.941;
          2, 45,0.948;
          2, 50,0.945;


          3,  5,0.946;
          3, 10,0.962;
          3, 15,0.944;
          3, 20,0.961;
          3, 25,0.953;
          3, 30,0.961;
          3, 35,0.953;
          3, 40,0.948;
          3, 45,0.936;
          3, 50,0.941;

          4,  5,0.957;
          4, 10,0.964;
          4, 15,0.960;
          4, 20,0.951;
          4, 25,0.956;
          4, 30,0.962;
          4, 35,0.956;
          4, 40,0.956;
          4, 45,0.933;
          4, 50,0.959;

          5,  5,0.928;
          5, 10,0.955;
          5, 15,0.948;
          5, 20,0.945;
          5, 25,0.950;
          5, 30,0.959;
          5, 35,0.965;
          5, 40,0.948;
          5, 45,0.953;
          5, 50,0.961;

          6,  5,0.943;
          6, 10,0.946;
          6, 15,0.946;
          %6, 20,0.915;
          6, 25,0.957;
          6, 30,0.954;
          6, 35,0.947;
          6, 40,0.962;
          6, 45,0.953;
          6, 50,0.963;


          %Ricordarsi che tutte e due le reti devono avere la stessa width

          
          ];


AUC_score = [        
           1,  5, 0.9429;
           1, 10, 0.9682;
           1, 15, 0.9688;
           1, 20, 0.9714;
           1, 25, 0.9601;
           1, 30, 0.9696;
           1, 35, 0.9654;
           1, 40, 0.9597;
           1, 45, 0.9721;
           1, 59, 0.9740;

           5,  5, 0.9815;
           5  10, 0.9799;
           5  15, 0.9937;
           5, 20, 0.9870;
           5, 25, 0.9956;
           5, 30, 0.9813;
           5, 35, 0.9959;
           5, 40, 0.9957;
           5, 50, 0.9892;

           10, 5, 0.9732;
           10,10, 0.9843;
           10,15, 0.9862;
           10,20, 0.9747;
           10,25, 0.9938;
           10,30, 0.9943;
           10,35, 0.9864;
           10,40, 0.9600;
           10,45, 0.9903;
           10,50, 0.9793;

           15, 5, 0.9808;
           %15,10, 0.9198;
           15,15, 0.9883;
           15,20, 0.9884;
           15,25, 0.9896;
           15,30, 0.9751;
           15,35, 0.9824;
           15,40, 0.9934;
           15,45, 0.9840;
           15,50, 0.9828;

           ];


%--------------------------------------------------------------------
%Accuracy
x=accuracy(:,1);
y=accuracy(:,2);
z=accuracy(:,3);
tri = delaunay(x,y);


figure(1)
trisurf(tri,x,y,z)
xlabel('N° hidden layer')
ylabel('Width')
zlabel('Accuracy')
title('Banchmark Accuracy')
colormap default
colorbar

figure(2)
trisurf(tri,x,y,z)
shading flat
xlabel('N° hidden layer')
ylabel('Width')
zlabel('Accuracy')
title('Banchmark Accuracy')
colormap jet
colorbar



figure(3)
trisurf(tri,x,y,z)
shading interp
xlabel('N° hidden layer')
ylabel('Width')
zlabel('Accuracy')
title('Banchmark Accuracy')
colormap jet
colorbar


%------------------------------------------------------------------------
%AUC_SCORE

x=AUC_score(:,1);
y=AUC_score(:,2);
z=AUC_score(:,3);
tri = delaunay(x,y);


figure(4)
trisurf(tri,x,y,z)
xlabel('N° hidden layer')
ylabel('Width')
zlabel('Accuracy')
title('Banchmark AUC SCORE')
colormap default
colorbar

figure(5)
trisurf(tri,x,y,z)
shading flat
xlabel('N° hidden layer')
ylabel('Width')
zlabel('Accuracy')
title('Banchmark AUC SCORE')
colormap jet
colorbar



figure(6)
trisurf(tri,x,y,z)
shading interp
xlabel('N° hidden layer')
ylabel('Width')
zlabel('Accuracy')
title('Banchmark AUC SCORE')
colormap jet
colorbar


