% Fig_point_alpha_lattice.m, 7.2.2025
% To show how points converge for point alpha

% Please cite as "Henning U. Voss and Douglas J. Ballon, Quasilattices of the aperiodic Spectre monotile, arXiv (2025)"
% The license attached in GitHub applies, at https://github.com/henningle/TileOneOne_Quasi

clear
close all

%% Parameters

write_figures=true;

Nmax=1;

figsize=500;

%% Tiling

[S,centers,xangles,vecs,N,Ncorners]=TileOneOne_fc(Nmax);

%% Paper figure

markersize_dataunits=1; % The x dimension of unit tile is about 4.4
markercolor=[222,105,54]/255;
C=cos(30*pi/180);

xys=[
    % -2.78, -1.48 % alpha from oservation, minimum of pairwise entropy
    % -2.7773, -1.4739 % alpha numetrically estimated below
    - (27*3^(1/2))/28 - 31/28, 3^(1/2)/28 - 43/28 % alpha symbolically estimated below
    ];

x=xys(1,1); y=xys(1,2);

h=figure('position',[100.,100.,figsize,figsize]);
ax=axes('Position', [0.1, 0.1, .8, .8]);

r=sqrt(x^2+y^2); phi=atan2(y,x);
points=[centers(:,1)+r*cos(phi+xangles),centers(:,2)+r*sin(phi+xangles)];

transparency=0;
plot(ax,S(:,1),S(:,2),'k',LineWidth=.2,Color=transparency*[1,1,1]);
hold on

tmpplot=plot(ax,points(:,1), points(:,2),'.','Marker','o','MarkerFaceColor',markercolor,'MarkerEdgeColor',markercolor,'MarkerSize',1);

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
axis off
axis image

% Plot connection arrow to center
for kk=1:N
    plot(ax,[centers(kk,1),points(kk,1)],[centers(kk,2),points(kk,2)],'k')
    plot(ax,centers(kk,1),centers(kk,2),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6)
    % text(ax,centers(kk,1)+.1,centers(kk,2),num2str(kk))
end
hold off

title('Point \alpha on first tiling iteration')

set(h,'Color', [1 1 1])

savefile=['fig_point_alpha'];
if write_figures
    print(h, '-dpng',  '-r400', [savefile '.png'])
    % saveas(h, [savefile '.fig'])
end

%% Analysis
% In-tile coordinate system like in paper

% disp(['x diameter = ' num2str(2*C+1.5) ', y diameter = ' num2str(0)])
% disp(['x midpoint = ' num2str(C+.75) ', y center = ' num2str(0)])
% disp(['x center = ' num2str(mean(-[0, 0, C, C-.5, C, 2*C, 2*C, 1+2*C, 2*C+1.5, C+1.5,C+1.5, C+.5, C, 0])) ...
%     ', y center = ' num2str(mean([0, 1, 1.5, 1.5+C, 1.5+2*C, 1+2*C, 2*C, 2*C, C, C-.5, C-1.5, C-1.5, -1.5, -1]))]) % Long form
disp(['x center = ' num2str(-(15*C+5.5)/14) ', y center = ' num2str((13*C+0.5)/14)]) % Short form

h=figure('position',[100.,100.,2.8*figsize,2.8*figsize]);
ax=axes('Position', [0.05, 0.05, .9, .9]); % Small margin

plot(ax,S(:,1),S(:,2),'k',LineWidth=1,Color=transparency*[1,1,1]);
hold on
tmpplot=plot(ax,points(:,1), points(:,2),'.','Marker','o','MarkerFaceColor',markercolor,'MarkerEdgeColor',markercolor,'MarkerSize',1);
set(tmpplot, 'MarkerSize', markersize_points);

% Plot connection arrow to center
for kk=1:N
    plot(ax,[centers(kk,1),points(kk,1)],[centers(kk,2),points(kk,2)],'k')
    plot(ax,centers(kk,1),centers(kk,2),'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6)
end

% Label vertices numerically
for kk=1:size(S,1) % 144
    if abs(S(kk,1))<1e-4; S(kk,1)=0; end
    if abs(S(kk,2))<1e-4; S(kk,2)=0; end
    text(S(kk,1)+.06,S(kk,2)-.06,['(' num2str(S(kk,1),3) ',' num2str(S(kk,2),3) ')'],'Color','b','FontSize',12)
end

% S1: Label vertices algebraically
xtexts={'0', '0', '-C', '-C+.5', '-C', '-2C', '-2C', '-2C-1', '-2C-1.5', '-C-1.5','-C-1.5', '-C-.5', '-C', '0'};
ytexts={'0', '1', '1.5', 'C+1.5', '2C+1.5', '2C+1', '2C', '2C', 'C', 'C-.5', 'C-1.5', 'C-1.5', '-1.5', '-1'};
for kk=1:14
    text(S(kk,1)+.06,S(kk,2)+.1,['(' xtexts{kk} ',' ytexts{kk} ')'],'FontSize',12)
end

% S3: Label vertices algebraically
xtexts={'-3C-3','-3C-3.5','-3C-3','-4C-3','-4C-3','-4C-2','-4C-1.5','-3C-1.5','-2C-1.5','-2C-2','-2C-1.5','-3C-1.5','-3C-1.5','-3C-2.5'};
ytexts={'1.5','-C+1.5','-2C+1.5','-2C+1','-2C','-2C', '-C','-C-.5','-C','0','C','C+.5','C+1.5','C+1.5'};
for kk=1:14
    kk2=2*16+kk;
    text(S(kk2,1)+.06,S(kk2,2)+.1,['(' xtexts{kk} ',' ytexts{kk} ')'],'FontSize',12)
end

% Test: Overlabel vertices numerically
%{
% First, put multipliers in
xtexts={'-3*C-3','-3*C-3.5','-3*C-3','-4*C-3','-4*C-3','-4*C-2','-4*C-1.5','-3*C-1.5','-2*C-1.5','-2*C-2','-2*C-1.5','-3*C-1.5','-3*C-1.5','-3*C-2.5'};
ytexts={'1.5','-C+1.5','-2*C+1.5','-2*C+1','-2*C','-2*C', '-C','-C-.5','-C','0','C','C+.5','C+1.5','C+1.5'};
for kk=1:14
       text(eval(xtexts{kk})+.06,eval(ytexts{kk})-.06,['(' num2str(eval(xtexts{kk}),3) ',' num2str(eval(ytexts{kk}),3) ')'],'Color','r','FontSize',12)
end
%}

% S4: Label vertices algebraically
xtexts={'-C-1.5','-C-2.5','-C-3','-2C-3','-3C-3','-3C-2.5','-3C-1.5','-3C-1.5','-2C-1.5','-2C-1','-2C','-2C','-C','-C-.5'};
ytexts={'3C+1.5','3C+1.5','2C+1.5','2C+2','2C+1.5','C+1.5', 'C+1.5','C+.5','C','2C','2C','2C+1','2C+1.5','3C+1.5'};
for kk=1:14
    kk2=3*16+kk;
    text(S(kk2,1)+.06,S(kk2,2)+.1,['(' xtexts{kk} ',' ytexts{kk} ')'],'FontSize',12)
end

% Test: Overlabel vertices numerically
%{
% First, put multipliers in
xtexts={'-C-1.5','-C-2.5','-C-3','-2*C-3','-3*C-3','-3*C-2.5','-3*C-1.5','-3*C-1.5','-2*C-1.5','-2*C-1','-2*C','-2*C','-C','-C-.5'};
ytexts={'3*C+1.5','3*C+1.5','2*C+1.5','2*C+2','2*C+1.5','C+1.5', 'C+1.5','C+.5','C','2*C','2*C','2*C+1','2*C+1.5','3*C+1.5'};
for kk=1:14
       text(eval(xtexts{kk})+.06,eval(ytexts{kk})-.06,['(' num2str(eval(xtexts{kk}),3) ',' num2str(eval(ytexts{kk}),3) ')'],'Color','r','FontSize',12)
end
%}

% Label centers numerically
% for kk=1:N
%     text(centers(kk,1)+.06,centers(kk,2),['(' num2str(centers(kk,1),3) ',' num2str(centers(kk,2),3) ')'],'Color','b')
%     text(ax,centers(kk,1)-.2,centers(kk,2),num2str(kk),'Color','r') % Tile number
% end

%% Label specific centers algebraically

% Note: angles=[150,90,90,30,-30,-30,-90] for N=1 means that tile 3 is rotated 150 and tile 4 90 degrees, because tile 2 is part of M
% [cos, -sin] [x]
% [sin,  cos] [y]
% Tile 3: 
% [cos(deg2rad(150)) -sin(deg2rad(150))] % = [-C -.5]
% [sin(deg2rad(150))  cos(deg2rad(150))] %   [.5 -C]
% Tile 4:
% [cos(deg2rad(90)) -sin(deg2rad(90))] % = [0 -1]
% [sin(deg2rad(90))  cos(deg2rad(90))] %   [1  0]

kk=1;
T10=[0;0]; % Origin of tile 1
C1=[-(15*C+5.5)/14; (13*C+0.5)/14]; % Center of tile 1
text(centers(kk,1)+.06,centers(kk,2)+.08, '-(15*C+5.5)/14','Color','r'); text(centers(kk,1)+.06,centers(kk,2)-0.08, '(13*C+0.5)/14','Color','r')

set(h,'Color', [1 1 1]) 
axis image
hold off

savefile='fig_point_alpha_shifts';
if write_figures
    print(h, '-dpng',  '-r400', [savefile '.png'])
    % saveas(h, [savefile '.fig'])
end

% kk=3;
% How to get center for tile 3
% The centers for tiles > 1 are not computed by the mean values of the vertices but by the shift and rotation of tile 1
% Use the relationship to the origin, rotate this vector by 90 degrees, and shift it by the new origin of tile 3
M3=[-C -.5; .5 -C]; % Rotation of tile 3
T30=[-3*C-3;1.5]; % Origin of tile 3
C3=M3*C1+T30; % Center of tile 3

% kk=4;
M4=[0 -1; 1 0];
T40=[-C-1.5;3*C+1.5];
C4=M4*C1+T40;

%% Determine point alpha as the common point for decoration of tiles 3 and 4

% Tile 1: alpha1=alpha relative to T01 = C1+[x;y]. Nothe this is in absolute coordinates in tile 1; should be [-4.1;-0.64] in absolute coords, where the decoration is
% Tile 2 is not needed (this is the lower body of the Mystic)
% Tile 3: alpha3=M3*(alpha relative to C1)+T30
% Tile 4: alpha4=M4*(alpha relative to C1)+T40

% First, test if present alpha fulfills this. The common point is [-1.72;0.0039]
% alpha3=M3*([x;y]+C1)+T30
% alpha4=M4*([x;y]+C1)+T40
% Yes, they are very similar but not exact

% Find exact value by setting alpha3=alpha4 and solving for x,y

zeta=inv(M4-M3)*(T30-T40);
% zeta=inv([C, -0.5; 0.5, C])*[-2*C-1.5;-3*C] % same
alpha_num=zeta-C1 %   -2.7773, -1.4739

% But what is the algebraic value?
% Do the exact same thing in symbolic algebra

Csym=sym(cos(pi/6)); % sqrt(3)/2 = cos30 = 0.8660

C1sym=[-(15*Csym+11/2)/14; (13*Csym+1/2)/14]; % Center of tile 1
M3sym=[-Csym -1/2; 1/2 -Csym]; % Rotation of tile 3
T30sym=[-3*Csym-3; 3/2]; % Origin of tile 3

M4sym=[0 -1; 1 0];
T40sym=[-Csym-3/2;3*Csym+3/2];

% Solving the inverse equation; resolve to alpha
zetasym=inv(M4sym-M3sym)*(T30sym-T40sym);
alpha_sym=zetasym-C1sym ;
alpha_sym=simplify(alpha_sym);
pretty(alpha_sym)
alpha_sym_eval=eval(alpha_sym)

% All the other merged points result from the same geometry modulo shift and rotation, as can be seen from the supplementary figure
% QED