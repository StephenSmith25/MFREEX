clear all
close all

PLOT_GRAPHS = true;
PLOT_DOMAINS_INFLUENCE = true; 


WITH_MOULD = false;
TMAX = 0.5;





path = './../../build/bin/preform/Displacement/';
pathSR = './../../build/bin/preform/srRod/';
boundaryNodes = csvread('./../../build/bin/preform/boundary.txt');
path_base = './../../build/bin/preform';

boundaryNodes = [boundaryNodes;boundaryNodes(1)];

addpath(path)
displacementdir = path ;
d = dir(displacementdir);
d1 = dir([displacementdir,'*.csv']);
numFiles = size(d,1) -3 ;

plotFiles = ceil(linspace(1,numFiles,10));

filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');
disp = csvread(filename,1);
boundaryNodes = boundaryNodes(5:80);


%figure
set(gcf, 'PaperPosition', [0 0 1 1]); %
set(gcf, 'PaperUnits', 'inches');
 x_width=7.25 ;y_width=9.125
 set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% -------------------------------------------------------------------------%
%                           PLOT 1
% -------------------------------------------------------------------------%


filename = strcat(path,'displacement_',num2str(plotFiles(1)),'.csv');


disp = csvread(filename,1);


hold on
rectangle('Position',[0 -150 100 280 ])
%hold on
%plot(0,110,'b*')
axis off
fill([disp(boundaryNodes,1);0],[disp(boundaryNodes,2);78.3800],'r','markersize',3)           % line plot
hold on
axis off
axis equal
plot([disp(boundaryNodes,1);0],[disp(boundaryNodes,2);78.3800],'g','markersize',3)           % line plot
xlim([0,32])
ylim([-220,85])


set(gcf, 'Position', [0 0 124 886]); %

export_fig './Barchart.png' -native
figure



A = imread('test.png');
B = imread('Barchart.png');

imshow(B);


size(B) % my image
size(A) 
RB = imref2d(size(B),0.43,0.43);

subplot(1,2,1)
imshow(B,RB)

RA = imref2d(size(A),0.26,0.26);

subplot(1,2,2)
imshow(A,RA)

figure
imshowpair(A, RA,B,RB, 'montage')


% -------------------------------------------------------------------------%
%                           PLOT 2
% -------------------------------------------------------------------------%


filename = strcat(path,'displacement_',num2str(plotFiles(5)),'.csv');


disp = csvread(filename,1);
figure

hold on
rectangle('Position',[0 -150 100 280 ])
%hold on
%plot(0,110,'b*')
axis off
fill([disp(boundaryNodes,1);0],[disp(boundaryNodes,2);78.3800],'r','markersize',3)           % line plot
hold on
axis off
axis equal
plot([disp(boundaryNodes,1);0],[disp(boundaryNodes,2);78.3800],'g','markersize',3)           % line plot
xlim([0,32])
ylim([-220,85])


set(gcf, 'Position', [0 0 124 886]); %

export_fig './Barchart.png' -native
figure



A = imread('test.png');
B = imread('Barchart.png');

imshow(B);


size(B) % my image
size(A) 
RB = imref2d(size(B),0.43,0.43);

subplot(1,2,1)
imshow(B,RB)

RA = imref2d(size(A),0.26,0.26);

subplot(1,2,2)
imshow(A,RA)

figure
imshowpair(A, RA,B,RB, 'montage')

