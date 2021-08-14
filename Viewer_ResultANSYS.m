%ANSYSviewer v1.0
%By Dew

clear all
clc

file = fopen('ansys cantilever.txt');
header = cell2mat( textscan(file,'%f',1,'collectoutput',1));
data = cell2mat( textscan(file, '%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter','\t','collectoutput',1));
fclose(file);

h = header;

% Plot Parameter
dispscale = 100;
scattersize = 20;
convertunits = 1;   % 1 to keep current units
                    % 2 to convert unit from mm to m


x = data(:,1);
y = data(:,2);
z = data(:,3);
vonmises = data(:,4);
sigmax = data(:,5);
sigmay = data(:,6);
sigmaxy = data(:,7);
Ux = data(:,8);
Uy = data(:,9);
Uz = data(:,10);

datasize = length(x);

for i = 1:datasize
    Xnew(i) = x(i) + Ux(i)*dispscale;
    Ynew(i) = y(i) + Uy(i)*dispscale;
    Znew(i) = z(i) + Uz(i)*dispscale;
end

if convertunits == 2
    for i = 1:datasize
        Xnew(i) = Xnew(i)/1000;
        Ynew(i) = Ynew(i)/1000;
        Znew(i) = Znew(i)/1000;
        sigmax(i) = sigmax(i)*1000;
        sigmay(i) = sigmay(i)*1000;
        sigmaxy(i) = sigmaxy(i)*1000;
        vonmises(i) = vonmises(i)*1000;
        Ux(i) = Ux(i)/1000;
        Uy(i) = Uy(i)/1000;
        Uz(i) = Uz(i)/1000;
        x(i) = x(i)/1000;
        y(i) = y(i)/1000;
        z(i) = z(i)/1000;
    end
end
    

figure()
set(gca,'FontSize', 18);
scatter3(Xnew,Ynew,Znew,scattersize,Uz,'filled')

colormap(jet);
daspect([1 1 1])
pbaspect([1 1 1])
h = rotate3d;
h.RotateStyle = 'box';
h.Enable = 'on';



% Probed Data Processing
file = fopen('geom_probe.txt');
probedata = cell2mat( textscan(file, '%u','Delimiter','\n','collectoutput',1));
fclose(file);

xprobe = [];
zprobe = [];

for i = 1:datasize
    if probedata(i) == 1
        xprobe(i) = x(i);
        zprobe(i) = Uz(i);
    end
end

[xprobe,sortIdx] = sort(xprobe,'ascend');
zprobe = zprobe(sortIdx);



