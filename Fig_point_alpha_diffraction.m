% Fig_point_alpha_diffraction.m, 7.2.2025
% Calculate the diffraction for Tile(1,1)

% Please cite as "Henning U. Voss and Douglas J. Ballon, Quasilattices of the Spectre monotile, arXiv (2025)"
% The license attached in GitHub applies, at https://github.com/henningle/TileOneOne_Quasi

clear
close all

%% Parameters

write_figures=true;

figsize=1000; % For spectral plot. < 1000 is a desaster
fontsize=12;
markercolor=[222,105,54]/255; % For plotting tiling

Nmax=5;

%% Tiling
[S,centers,xangles,vecs,N,Ncorners]=TileOneOne_fc(Nmax);

%% Main

markersize_dataunits=2.6; % The x dimension of unit tile is about 4.4

xys=[
    -3 1.8
    ];

x=xys(1,1); y=xys(1,2);
r=sqrt(x^2+y^2); phi=atan2(y,x);

%% Masked pattern

% Circular ROI to avoid artifacts
switch Nmax
    case 4
        fc=[6,44]; % Field center for Nmax=4
        radius=44;
        scattersize=120;
    case 5
        fc=[4,128]; % Field center for Nmax=4
        radius=130;
        scattersize=120/3;
end
centers2=centers; % The centers that will actually be used for spectrum
for ii=N:-1:1
     if norm(centers2(ii,:)-fc)>radius; centers2(ii,:)=[]; xangles(ii)=[]; end
end
N2=size(centers2,1);
disp(['Cells for Tile(1,1)=' num2str(N2/7)])
points2=[centers2(:,1)+r*cos(phi+xangles),centers2(:,2)+r*sin(phi+xangles)]; % xangles is vector for all centers

figure('position',[400.,100.,figsize,figsize]);
ax=axes('Position', [0, 0, 1, 1]); % Fill the whole frame, no margin, no border. Note next line, too

% Initial plot with arbitrary MarkerSize
tmpplot=plot(ax,points2(:,1), points2(:,2),'.','Marker','o','MarkerFaceColor',markercolor,'MarkerEdgeColor',markercolor,'MarkerSize',1);
hold on

% Correct scaling of markers. This is needed for defining it in terms of Spectre units
x_range=diff(xlim(gca)); y_range=diff(ylim(gca));
fig_position=get(gcf, 'Position'); % in pixels [left bottom width height]
ax_position=get(gca, 'Position');  % normalized [left bottom width height]
fig_width_inch=fig_position(3) * ax_position(3) / get(gcf, 'ScreenPixelsPerInch');
fig_height_inch=fig_position(4) * ax_position(4) / get(gcf, 'ScreenPixelsPerInch');
x_marker_size_points=markersize_dataunits / x_range * fig_width_inch * 72;
y_marker_size_points=markersize_dataunits / y_range * fig_height_inch * 72;
markersize_points=mean([x_marker_size_points, y_marker_size_points]);
set(tmpplot, 'MarkerSize', markersize_points);
axis image

%% Magnitude of the Fourier transform

% Continuous FT
kxs=linspace(-pi, pi, figsize); 
kys=linspace(-pi, pi, figsize);
[kx_grid, ky_grid]=meshgrid(kxs, kys);

% Sum contributions of each delta peak with index ii with tapering amplitude A(ii)
A=ones(N2,1); % No tapering
F=zeros(size(kx_grid));
for ii=1:N2
    F=F+A(ii)*exp(-1j*(kx_grid*points2(ii,1)+ky_grid*points2(ii,2)));
end
F=abs(F); F=F/max(F(:));

% lambda=2*pi/norm([(-0.0408816+0.00314474),-0.789329-0.00314474])

figure('position',[400.,100.,figsize,figsize]);
ax=axes('Position', [0, 0, 1, 1]);
imagesc(ax,kxs, kys, F.^.5); % As always, use ' for imagesc
axis xy
axis square
axis off

colormap jet
if write_figures
    print(gcf, '-dpng',  '-r400', [ 'png/fig_point_alpha_diffraction_Nmax' num2str(Nmax) '.png'])
end
