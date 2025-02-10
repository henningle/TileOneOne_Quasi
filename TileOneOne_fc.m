% TileOneOne_fc.m, 5.12.2024

% Please cite as "Henning U. Voss and Douglas J. Ballon, Quasilattices of the aperiodic Spectre monotile, arXiv (2025)"
% The license attached in GitHub applies, at https://github.com/henningle/TileOneOne_Quasi

function [S,centers,xangles,vecs,N,Ncorners]=TileOneOne_fc(Nmax)

%% Parameters

plotting=false; % Plotting quasilattice
plotkeypoints=false; % Plotting key points

%% Main script

% Tile placement function
% Each tile has 3 components: x, y, key point number (0 for no key point, 1 to 4 if key point)
% The placetile function takes in a matrix of size [number of rows, 3]
% The first 2 columns are the coordinates of the tiling, the third is a label of key points
placetile=@(x, angle, coord) ...
    [transpose( [cos(angle*pi/180), -sin(angle*pi/180); ...
    sin(angle*pi/180), cos(angle*pi/180)]*x(:,1:2)' ...
    +transpose(repmat(coord,[size(x,1),1]))) , x(:,3)];

mirror=@(x) [transpose([-1, 0; 0, 1]*(x(:,1:2)')),  x(:,3)];

Ncorners=15; % 14 vertices + 1 vertex to get back to origin

% Number of expected tiles for each iteration
nS=zeros(1,Nmax+1);
nM=nS;
nS(1)=1; nM(1)=2; % This is index 0
for n=2:Nmax+1
    nS(n)=nM(n-1)+7*nS(n-1);
    nM(n)=nM(n-1)+6*nS(n-1);
end
nS=nS(2:end);
% nM=nM(2:end);
disp(['Predicted number of tiles for N = 1, 2, ...: ' num2str(nS)])

%% S0 and M0

% Some definitions
cos30=cos(30*pi/180); sin30=sin(30*pi/180);
cos60=cos(60*pi/180); sin60=sin(60*pi/180);
cos120=cos(120*pi/180); sin120=sin(120*pi/180);
cos150=cos(150*pi/180); sin150=sin(150*pi/180);
cos210=cos(210*pi/180); sin210=sin(210*pi/180);
cos240=cos(240*pi/180); sin240=sin(240*pi/180);
cos300=cos(300*pi/180); sin300=sin(300*pi/180);
cos330=cos(330*pi/180); sin330=sin(330*pi/180);

% One Tile(1,1)
% Each tile has 14 edges and 15 points, the last point for closure
% Plus one NaN for separation of tiles while plotting
% x-components of specter S
Stmp1=[ 0
    0
    cos30
    cos30+cos120
    cos30+cos120+cos60
    cos30+cos120+cos60+cos330
    cos30+cos120+cos60+cos330
    cos30+cos120+cos60+cos330+1
    cos30+cos120+cos60+cos330+1+cos300
    cos30+cos120+cos60+cos330+1+cos300+cos210
    cos30+cos120+cos60+cos330+1+cos300+cos210
    cos30+cos120+cos60+cos330+cos300+cos210
    cos30+cos120+cos60+cos330+cos300+cos210+cos240
    cos30+cos120+cos60+cos330+cos300+cos210+cos240+cos150
    cos30+cos120+cos60+cos330+cos300+cos210+cos240+cos150
    NaN];
% y-components of specter S
Stmp2=[0
    1
    1+sin30
    1+sin30+sin120
    1+sin30+sin120+sin60
    1+sin30+sin120+sin60+sin330
    sin30+sin120+sin60+sin330
    sin30+sin120+sin60+sin330
    sin30+sin120+sin60+sin330+sin300
    sin30+sin120+sin60+sin330+sin300+sin210
    sin30+sin120+sin60+sin330+sin300+sin210-1
    sin30+sin120+sin60+sin330+sin300+sin210-1
    sin30+sin120+sin60+sin330+sin300+sin210-1+sin240
    sin30+sin120+sin60+sin330+sin300+sin210-1+sin240+sin150
    sin30+sin120+sin60+sin330+sin300+sin210+sin240+sin150
    NaN];
% Specter S
S=[Stmp1,Stmp2];

% Augment S with key points to be carried over during iterative construction
keypoints=zeros(16,1); keypoints(4)=1; keypoints(6)=2; keypoints(8)=3; keypoints(14)=4;
S=[S,keypoints];
keyind=find(S(:,3)>0)'; % The indexes that point to key points
% S(keyind,3)' % Should be 1,2,3,4 for iteration 0

% The Mystic
M=[S; placetile(S,-30,[sqrt(3)/2,-1.5])];
M(:,3)=0; % M should not have key points not to contaminate S

%% Iteration loop

for N=1:Nmax

    S=mirror(S);
    M=mirror(M);

    % The first two iterations require explicit key points for the next iteration
    switch N
        case 1
            Sind=[4,6,8,14]; % Key points on S0, for making S1
            Mind=Sind+16;
            angles=[150,90,90,30,-30,-30,-90];
        case 2
            Sind=[84,118,132,38]; % Key points on S1, for making S2
            % Indexing changes in M2 due to the missing patch
            Mind=[84-16,118-16,132-16,38]; %
            angles=[120,60,60,0,-60,-60,-120]; % Conjecture: Angles converged from here on
    end

    % For iteration 3 and higher, key points are computed from the previous two iterations
    if N > 2
        Sind=circshift(keyind,-1);
        Mind(1:3)=Sind(1:3)-16*nS(N-2);
    end

    % Define some tiling cells
    cell0=M;
    cell1=placetile([S(:,1:2)-S(Sind(3),1:2), S(:,3)],angles(1), M(Mind(3),1:2));
    cell2=placetile([S(:,1:2)-S(Sind(2),1:2), S(:,3)],angles(2),cell1(Sind(4),1:2));
    cell3=placetile([S(:,1:2)-S(Sind(3),1:2), S(:,3)],angles(3),cell2(Sind(1),1:2));
    cell4=placetile([S(:,1:2)-S(Sind(2),1:2), S(:,3)],angles(4),cell3(Sind(4),1:2));
    cell5=placetile([S(:,1:2)-S(Sind(2),1:2), S(:,3)],angles(5),cell4(Sind(4),1:2));
    cell6=placetile([S(:,1:2)-S(Sind(3),1:2), S(:,3)],angles(6),cell5(Sind(1),1:2));
    cell7=placetile([S(:,1:2)-S(Sind(2),1:2), S(:,3)],angles(7),cell6(Sind(4),1:2));

    %% Assemble S and M

    S=[cell0 % 2x16 elements because this is the mystic M
        cell1 % 16 elements
        cell2
        cell3
        cell4
        cell5
        cell6
        cell7
        ];
    M=[cell0 % 2x16 elements
        cell1 % 16 elements
        cell2
        cell4 % Note we skipped cell3
        cell5
        cell6
        cell7
        ];
    M(:,3)=0; % Again, M should not have key points

    if plotting
        figure('position',[500.,100.,1000,1000]);
        plot(S(:,1),S(:,2),'k',LineWidth=.1);
        axis image
        title(['S' num2str(N)])
    end

    %% Key points of this iteration

    keyind=find(S(:,3)>0); % Always 7*4=28 elements after pruning previous iteration

    if plotting
        if plotkeypoints
            % Plot all key points before pruning, in red
            hold on
            for n=1:4:length(keyind)
                plot(S([keyind(n:n+3);keyind(n)],1),S([keyind(n:n+3);keyind(n)],2),'ro-');
            end
        end
    end

    %% Key points for next iteration

    % At the beginning, there are always 9*4=36 key points. At the end, there should be 4
    keyind2=keyind([2,13,22,25]); % Read off from Smith23b Fig. A1.2a and previous plot
    if N>1; keyind2=keyind([2,13,22,25]+1); end % Don't ask
    % S(keyind2,3)' % Should be (2,1,2,1) for iteration 1 and (3,2,3,2) for iteration >= 2

    S(:,3)=0; % Clean out this iteration's key points

    % Set new key points for constructing the next iteration
    S(keyind2(1),3)=1;
    S(keyind2(2),3)=2;
    S(keyind2(3),3)=3;
    S(keyind2(4),3)=4;

    keyind=find(S(:,3)>0); % Must be 38 84 118 132 for iteration 2

    if plotting
        if plotkeypoints
            % Plot only key points taken over to next iteration, in blue
            plot(S([keyind;keyind(1)],1),S([keyind;keyind(1)],2),'bo-');
        end
    end

    % Count tiles
    nStrue=size(S,1)/16;
    disp(['Iteration ' num2str(N) ' has ' num2str(nStrue) ' tiles'])

end

% Make the turtle always look to the left. It has to be inverted for even Nmax
if not(mod(Nmax,2))
    disp('Flipping pattern for even Nmax')
    S=mirror(S);
end

S=S(:,1:2); % Saving only the lattice proper
N=nStrue; % Size or RF array

%% Further outputs

centers=zeros(N,2);
for n=1:N % Loop through resonators
    centers(n,:)=[mean(S((n-1)*16+1:(n-1)*16+14,1)),mean(S((n-1)*16+1:(n-1)*16+14,2))]; % When computing the center, the last point needs to be removed becasue it is repeated
end

% Angles
% Paralell to x axis to make the angles multiples of 30 degrees
% 2 vectors pointing to the base edge define the base vector, which serves as a reference for the angle
% Astonishingly, there is no need to correct for flipped S for odd Nmax
vecs=NaN*zeros(3*N,2); % Has NaNs for easier plotting
xangles=zeros(N,1);
for ii=1:N
    vecs(3*ii-2,:)=[S((ii-1)*(Ncorners+1)+11,1),S((ii-1)*(Ncorners+1)+11,2)]; % vector 1
    vecs(3*ii-1,:)=[S((ii-1)*(Ncorners+1)+12,1),S((ii-1)*(Ncorners+1)+12,2)]; % vector 2
    % Angles
    v1=vecs(3*ii-2,:); v2=vecs(3*ii-1,:);
    vx=v2(1)-v1(1); vy=v2(2)-v1(2); % First get direction vector from v1 to v2
    xangles(ii)=atan2(vy, vx);
end
% rad2deg(xangles(1:10))' % Must start with 0