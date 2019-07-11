
close all
clear all




nnx = 31;

nny = 21;

h = 0.1;
w = 1;

x = linspace(0,w,nnx);
y = linspace(0,h,nny);

[X,Y] = meshgrid(x,y);

nodes = [X(:),Y(:)];
figure
plot(nodes(:,1),nodes(:,2),'bo');
K = convhull(nodes(:,1),nodes(:,2));



axis equal


for i = 1:length(K)-1
   segments(i,:) = [K(i),K(i+1),2] ;
    
end

for i = 1:length(segments)
   
    hold on
    plot(nodes(segments(i,1:2),1),nodes(segments(i,1:2),2),'r-');
    
end


hold on






%% write files
dlmwrite('../../problems/rubber_comp/rubber.nodes',[length(nodes),0],'delimiter',' ')
dlmwrite('../../problems/rubber_comp/rubber.nodes',nodes,'-append','delimiter',' ')


%segments
dlmwrite('../../problems/rubber_comp/rubber.segs',length(segments),'delimiter',' ')
dlmwrite('../../problems/rubber_comp/rubber.segs',segments,'-append',.........
    'delimiter',' ')

