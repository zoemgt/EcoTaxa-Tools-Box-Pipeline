
%% ── 12. Multivariate analyses (optional) ──────────────────────────────────
% Operates on the instrument-agnostic base(i).regroupped / data_ab_* /
% data_bv_* / plankton_groups_data_* structures built in EcoTaxa_Tool_Box.m,
% so it works unchanged for IFCB, FlowCam, PlanktoScope and ZooScan alike.

% Colormap used for dendrograms and NBSS-slope overlays (load once, reuse
% if already loaded by the NBSS section of EcoTaxa_Tool_Box.m)
if ~exist('slandarerCM', 'var')
    load('slanCM_Data.mat');
end

%% ── 12.1 Choose data and normalisation ─────────────────────────────────────

%Choose data to analyse
list3_a = {'Abundance per taxo','Abundance per plankton group','Abundance per trophic group','Biovolume per taxo','Biovolume per plankton group','Biovolume per trophic group'}; 
[indx3,tf] = listdlg('PromptString','What data should be used for analysis?','ListString',list3_a,'SelectionMode','single','ListSize',[400,200]);
list3 = {'data_ab_taxo','data_ab_group','data_ab_troph','data_bv_taxo','data_bv_group','data_bv_troph'}; 

    dataACP_name = list3{indx3};
    data = eval(dataACP_name);

    if indx3 == 1 || indx3 == 4
    group=plankton_groups_data_taxo;
    elseif indx3 == 2 || indx3 == 5
    group=plankton_groups_data_group;
    elseif indx3 == 3 || indx3 == 6
    group=plankton_groups_data_troph_name;
    end

    for i = 1:numel(group)
    group{i} = strrep(group{i}, '_', ' ');
    end

%Choose normalisation
list3 = {'No normalisation',  'Log x+1, more weight to rare', 'Double cube root, more weight to rare', 'Relative abundance', 'Hellinger transform, relative + more weight to rare', 'Centred and reduced, forced normalisation loss of relative strength', 'Box cox, crude forced normalisation loss of relative strength'}; 
[normalisation,tf] = listdlg('PromptString','With what normalisation?','ListString',list3,'SelectionMode','single','ListSize',[600,200]);

if normalisation==1; %No
data_stat=data';

elseif normalisation==2; %logx+1
data_stat=log(data+1)';

elseif normalisation==3;
data_stat=((data').^0.25);

elseif normalisation==4; %relative abundance
total=sum(data,2);
[n,m]=size(data);
total=total*ones(1,m);
relativeabund=100*data./total;
data_stat=relativeabund';

elseif normalisation==5; %hellinger transform
total=sum(data,2);
[n,m]=size(data);
total=total*ones(1,m);
relativeabund=100*data./total;
data_stat=(relativeabund').^0.5;

elseif normalisation==6;
data_stat=normalize(data);
data_stat=data_stat';

elseif normalisation==7;
for i=1:size(data,2);
[transdat,lambda] = boxcox(data(:,i)+0.0000001);
transdattemp(:,i)=transdat;
end
data_stat=normalize(transdattemp);
data_stat=data_stat';

end

%to save
parentDir = fullfile(pwd, 'multivariates analysis');
mkdir(parentDir);

%% ── 12.2 PCA ──────────────────────────────────────────────────────────────

%to save
subDir = 'PCA';
fullSubDirPath = fullfile(parentDir, subDir);
    mkdir(fullSubDirPath);

[pc,zscores,latent,tsquared,pcvars] = pca(data_stat','Algorithm','svd'); 

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
zscores = (zscores./repmat(max(max(abs(zscores))), size(zscores)));

% PCA plot: 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

scatter(zscores(:,1),zscores(:,2),50,'filled','MarkerFaceColor', [0.9 0.9 0.9],'markeredgecolor',[0.5 0.5 0.5])

% Plot in function of the number of groups 

if numel(group) <= 15

[p,q]=size(pc);
hold on
for i=1:p
    if ((pc(i,1)^2+pc(i,2)^2)^0.5)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

elseif numel(group) > 15

[p,q]=size(pc);
hold on

distances = sqrt(pc(:,1).^2 + pc(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

end

xlabel(['PCA 1 ' num2str(fix(pcvars(1))) '% variance'],'fontsize',10);
ylabel(['PCA 2 ' num2str(fix(pcvars(2))) '% variance'],'fontsize',10);

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCA on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

% Representation with RGB colors 
% Extract RGB color on the PCA

    RGB=zscores(:,1:3);
    minRGB=min(RGB);
    RGB=RGB-ones(size(data_stat,2),1)*minRGB;
    maxRGB=max(RGB);
    RGB=RGB./(ones(size(data_stat,2),1)*maxRGB);
    RGB=1-RGB;

    axis1=1;
    axis2=2;
    axis3=3;

% Plot PCA with RGB colors 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

scatter(zscores(:,axis1),zscores(:,axis2),50,RGB,'filled','markeredgecolor',[0.5 0.5 0.5])

[p,q]=size(pc);
pc2=pc;
hold on

if numel(group) <= 15

[p,q]=size(pc);
hold on
for i=1:p
    if ((pc(i,1)^2+pc(i,2)^2)^0.5)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

elseif numel(group) > 15

[p,q]=size(pc);
hold on

distances = sqrt(pc(:,1).^2 + pc(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

end

xlabel(['PCA 1 ' num2str(fix(pcvars(1))) '% variance'],'fontsize',10);
ylabel(['PCA 2 ' num2str(fix(pcvars(2))) '% variance'],'fontsize',10);

% Map with RGB colors 

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);

m_scatter(lon,lat, 50, RGB, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCA and RGB color on map ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

%% ── 12.3 PCoA ─────────────────────────────────────────────────────────────

%to save
subDir = 'PCoA';
fullSubDirPath = fullfile(parentDir, subDir);
    mkdir(fullSubDirPath);

%Methods for matrice (lots of other choices in f_dis Fathom)
list6 = {'Euclidean, recommanded for Hellinger normalisation','Bray curtis, recommanded for others normalisation'}; 
[indx6,tf] = listdlg('PromptString','Distance methods for PCoA?','ListString',list6,'SelectionMode','single','ListSize',[400,200]);

if indx6==1
dis = f_dis(data_stat','euc');  
elseif indx6==2
dis = f_dis(data_stat','bc');  
end

pcoa = f_pcoa_modified(dis,1); 

% Scale scores to fit within -1 to 1 bounds (after eq. 1.10 in L&L, 1998):
scores = pcoa.scores;
expl   = pcoa.expl;
scores = (scores./repmat(max([max(abs(scores(:,1)));max(abs(scores(:,2)))]), size(scores)));

%PCOA plot (adapted from f_pcoaPlot (Fathom): 

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

scatter(scores(:,1), scores(:,2), 50, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'markeredgecolor', [0.5 0.5 0.5]);
%text(scores(:,1),scores(:,2),sampleid,'HorizontalAlignment','center', 'VerticalAlignment','middle','Color','k','Interpreter','tex'); 

Y=data_stat';
ncY = size(Y,2);
s = size(scores,2);
vec = zeros(ncY,s);
   for i = 1:s
      for j = 1:ncY
       vec(j,i) = f_corr(Y(:,j),scores(:,i))*std(Y(:,j))/std(scores(:,i));
      end
   end
vec = (vec./repmat(max([max(abs(vec(:,1)));max(abs(vec(:,2)))]), size(vec)));

eDis = f_dis([0 0;vec(:,1:2)],'euc');
eDis = eDis(2:end,1);

% Set delta inversely proportional to length of longest vector:
ratio = (abs(eDis./ repmat(max(eDis),size(eDis))));
ratio = ones(size(ratio))./ratio;
delta = 1 + (1* repmat(0,size(eDis)) .* ratio);
    
hold on 

if numel(group) <= 15

       n = size(vec,1); % get # of variables
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
   end

elseif numel(group) > 15

           n = size(vec,1); % get # of variables

distances = sqrt(vec(:,1).^2 + vec(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.4 * numel(distances)));

for j=topIndices'
    if distances(j)
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
    end
end

end

axis1 = sprintf('%2.2f',expl(1,1));
axis2 = sprintf('%2.2f',expl(2,1));
xlabel(['PCoA 1 (' num2str(axis1) ' % variance)'],'fontsize',10);
ylabel(['PCoA 2 (' num2str(axis2) ' % variance)'],'fontsize',10);

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCoA on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% PCoA plot representation with RGB

%Representation with RGB colors 
%Extract RGB color on the PCoA

    RGB=scores(:,1:3);
    minRGB=min(RGB);
    RGB=RGB-ones(size(data_stat,2),1)*minRGB;
    maxRGB=max(RGB);
    RGB=RGB./(ones(size(data_stat,2),1)*maxRGB);
    RGB=1-RGB;

    axis1=1;
    axis2=2;
    axis3=3;

%Plot PCoA with RGB colors

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

scatter(scores(:,1), scores(:,2), 50, RGB, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
%text(scores(:,1),scores(:,2),sampleid,'HorizontalAlignment','center', 'VerticalAlignment','middle','Color','k','Interpreter','tex'); 

Y=data_stat';
ncY = size(Y,2);
s = size(scores,2);
vec = zeros(ncY,s);
   for i = 1:s
      for j = 1:ncY
       vec(j,i) = f_corr(Y(:,j),scores(:,i))*std(Y(:,j))/std(scores(:,i));
      end
   end
vec = (vec./repmat(max([max(abs(vec(:,1)));max(abs(vec(:,2)))]), size(vec)));

eDis = f_dis([0 0;vec(:,1:2)],'euc');
eDis = eDis(2:end,1);

% Set delta inversely proportional to length of longest vector:
ratio = (abs(eDis./ repmat(max(eDis),size(eDis))));
ratio = ones(size(ratio))./ratio;
delta = 1 + (1* repmat(0,size(eDis)) .* ratio);
    
hold on 

if numel(group) <= 15

       n = size(vec,1); % get # of variables
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
   end

elseif numel(group) > 15

           n = size(vec,1); % get # of variables

distances = sqrt(vec(:,1).^2 + vec(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for j=topIndices'
    if distances(j)
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
    end
end

end

axis1 = sprintf('%2.2f',expl(1,1));
axis2 = sprintf('%2.2f',expl(2,1));
xlabel(['PCoA 1 (' num2str(axis1) ' % variance)'],'fontsize',10);
ylabel(['PCoA 2 (' num2str(axis2) ' % variance)'],'fontsize',10);

%Map with RGB colors

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);

m_scatter(lon,lat, 50, RGB, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

%Export and save? 
if doSave==false;
elseif doSave==true;
        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCoA and RGB color on map ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

%% ── 12.4 Correspondence analysis (CA) ────────────────────────────────────

%to save
subDir = 'CA';
fullSubDirPath = fullfile(parentDir, subDir);
    mkdir(fullSubDirPath);

% CA (inspired from Rencher, A. C. (2000))

%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A., R. Hernandez-Walls and K. Barba-Rojo. (2008). corran:
%    Correspondence Analysis. A MATLAB file. [WWW document]. URL http://
%    www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=21505
%  Reference:
%  Rencher, A. C. (2000), Methods of Multivariate Analysis. 2nd ed. 
%       John Wiley & Sons. p. 514-525

[a, b] = size(data_stat); 
N = sum(data_stat(:)); % total of the matrix
P = data_stat / N; % Correspondence matrix (relative frequencies)
r = sum(P, 2); % Row totals (marginal frequencies)
c = sum(P, 1); % Column totals (marginal frequencies)
Dr = diag(r); % Diagonal matrix of row totals
Dc = diag(c); % Diagonal matrix of column totals
R = inv(Dr) * P; % Row profiles (standardized)
C = P * inv(Dc); % Column profiles (standardized)
Z = sqrt(inv(Dr)) * (P - r * c) * sqrt(inv(Dc)); % Matrix for SVD
[U, S, V] = svd(Z); % Singular Value Decomposition
k = rank(inv(Dr) * (P - r * c) * inv(Dc) * (P - r * c)'); % Maximum number of dimensions
SV = diag(S); %singular values
SV = SV(1:k); %singular value equal to zero is deleted
EV = SV.^2; %eigenvalues (inertia)
PI = EV./sum(EV)*100; %percent of inertia
PA = cumsum(PI); %cumulative percent of inertia
x2 = N*EV; %chi-squared
A = sqrt(Dr)*U;
B = sqrt(Dc)*V;
X2 = N*trace(inv(Dr)*(P-r*c)*inv(Dc)*(P-r*c)'); %total chi-squared
v = (a-1)*(b-1); %degrees of freedom
P = 1-chi2cdf(X2,v); %P-value

% Scores for rows and columns
X=inv(Dr)*A*S; %row coordinates
X=X(:,1:k);
Y=inv(Dc)*B*S'; %column coordinates
Y=Y(:,1:k);

% CA plot:

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

% distances
distances_X = sqrt(X(:,1).^2 + X(:,2).^2);
distances_Y = sqrt(Y(:,1).^2 + Y(:,2).^2);

if numel(group) <= 15

    %scatter(X(:,1), X(:,2), 50, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(:,1), X(:,2), group', 'Color', 'k', 'fontsize', 8);

    hold on
    for i=1:size(X,1)
    plot([0 X(i,1)],[0 X(i,2)],'k-')
    end

elseif numel(group) > 15

    [~, sortedIndices_X] = sort(distances_X, 'descend');
    topIndices_X = sortedIndices_X(1:round(0.3 * numel(distances_X)));
    
    %scatter(X(topIndices_X, 1), X(topIndices_X, 2), 50, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(topIndices_X, 1), X(topIndices_X, 2), group(topIndices_X), 'Color', 'k', 'fontsize', 8);

    hold on
    for i = topIndices_X'
        plot([0 X(i,1)], [0 X(i,2)], 'k-');
    end

end

hold on;

scatter(Y(:,1), Y(:,2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);

xlim([-max(abs([X(:,1); Y(:,1)])) * 1.1, max(abs([X(:,1); Y(:,1)])) * 1.1]);
ylim([-max(abs([X(:,2); Y(:,2)])) * 1.1, max(abs([X(:,2); Y(:,2)])) * 1.1]);

%text(Y(:,1) + .01, Y(:,2) + .1, sampleid', 'Color', 'k', 'fontsize', 8);

xlabel(['CA 1 (' num2str(fix(PI(1))) '% variance)']);
ylabel(['CA 2 (' num2str(fix(PI(2))) '% variance)']);

%Export and save? 
if doSave==false;
elseif doSave==true;
    outputDir = fullfile(parentDir, subDir);
    mkdir(outputDir);
    
    filename = fullfile(outputDir, sprintf(['CA on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1);
close(gcf);

% Plot CA with RGB colors

%Representation with RGB colors 
%Extract RGB color on the CA

    RGB=Y(:,1:3);
    minRGB=min(RGB);
    RGB=RGB-ones(size(data_stat,2),1)*minRGB;
    maxRGB=max(RGB);
    RGB=RGB./(ones(size(data_stat,2),1)*maxRGB);
    RGB=1-RGB;

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

% distances
distances_X = sqrt(X(:,1).^2 + X(:,2).^2);
distances_Y = sqrt(Y(:,1).^2 + Y(:,2).^2);

if numel(group) <= 15

    %scatter(X(:,1), X(:,2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(:,1), X(:,2), group', 'Color', 'k', 'fontsize', 8);

    hold on
    for i=1:size(X,1)
    plot([0 X(i,1)],[0 X(i,2)],'k-')
    end

else

    [~, sortedIndices_X] = sort(distances_X, 'descend');
    topIndices_X = sortedIndices_X(1:round(0.3 * numel(distances_X)));
    
    %scatter(X(topIndices_X, 1), X(topIndices_X, 2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(topIndices_X, 1), X(topIndices_X, 2), group(topIndices_X), 'Color', 'k', 'fontsize', 8);

    hold on
    for i = topIndices_X'
        plot([0 X(i,1)], [0 X(i,2)], 'k-');
    end

end

hold on;

scatter(Y(:,1), Y(:,2), 50, RGB, 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5]);

xlim([-max(abs([X(:,1); Y(:,1)])) * 1.1, max(abs([X(:,1); Y(:,1)])) * 1.1]);
ylim([-max(abs([X(:,2); Y(:,2)])) * 1.1, max(abs([X(:,2); Y(:,2)])) * 1.1]);

%text(Y(:,1) + .01, Y(:,2) + .1, sampleid', 'Color', 'k', 'fontsize', 8);

xlabel(['CA 1 (' num2str(fix(PI(1))) '% variance)']);
ylabel(['CA 2 (' num2str(fix(PI(2))) '% variance)']);

%Map with RGB colors

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);

m_scatter(lon,lat, 50, RGB, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

%Export and save? 
if doSave==false;
elseif doSave==true;
    outputDir = fullfile(parentDir, subDir);
    mkdir(outputDir);
    
    filename = fullfile(outputDir, sprintf(['CA and RGB color on map on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1);
close(gcf);

%% ── 12.5 Hierarchical clustering ─────────────────────────────────────────

%to save
subDir = 'Clusters analysis';
fullSubDirPath = fullfile(parentDir, subDir);
    mkdir(fullSubDirPath);

% Dendrogram to define cut

% Dendrocut value suggestion
Z = linkage(data_stat', 'complete', 'euclidean'); 
incon=inconsistent(Z);
test=cumsum(flipud(incon(:,4)));
I=find(test==max(test));
dendrocut_auto=incon(I(1),1); %Suggestion

% Plot dendrogram with dendrocut_auto

if ~iscell(sampleid); sampleid = num2cell(sampleid); end
if ischar(sampleid) || isstring(sampleid); sampleid = cellstr(sampleid); end

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot(1,2,1)
plot(incon(:,1),cumsum(flipud(incon(:,4))),'Color', [0.5 0.5 0.5])
hold on
plot(incon(I(1),1),cumsum(flipud(incon(I(1),4))),'o','MarkerFaceColor', [0.93,0.69,0.13])
text(incon(I(1),1), cumsum(flipud(incon(I(1),4)))+5, sprintf('%.2f', dendrocut_auto), 'Color', [0.93,0.69,0.13], 'HorizontalAlignment', 'left');

dis2=pdist(data_stat', 'euclidean');
leafOrder = optimalleaforder(Z,dis2);

subplot(1,2,2)
[lineH,T, Perm] = dendrogram(Z,0, 'ColorThreshold',dendrocut_auto,'reorder',leafOrder,'Labels',sampleid);
colors = cell2mat(get(lineH,'Color'));
set(lineH,'Color', [0.5 0.5 0.5]);
set(gca, 'XTick', []);

hold on;
hline = refline(0,dendrocut_auto);
hline.Color = [0.93,0.69,0.13]; 

color = slandarerCM(6).Colors{1,9};  % spectral (change manually)

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir, subDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['Dendrogram for clustering on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

% Or set the dendrocut value manually?
list2 = {'Keep', 'Change'};  
[dendro,tf] = listdlg('PromptString','Keep automatic dendrocut or select manually?','ListString',list2,'SelectionMode','single','ListSize',[400,200]);

if dendro==1;
    dendrocut = dendrocut_auto; 
elseif dendro==2;
        name = 'Dendrocut value';
        prompt = {'Enter the Dendrocut value:'};
        dims = [1 35];
        definput = {'10'};
        answerfrac1 = inputdlg(prompt, name, dims, definput);
        dendrocut = str2double(answerfrac1{1});
end

close gcf

% Create and Plot Clustergram

CGobj=clustergram(data_stat,'RowLabels',group','ColumnLabels',sampleid,'Linkage','complete','ColumnPDist', 'euclidean','ColumnPDist', 'euclidean','Colormap',flipud(color),'Standardize',3,'Symmetric','false','DisplayRatio',[0.2 0.000000000000001],'OptimalLeafOrder',1);
CGobj.Dendrogram = dendrocut;
CGobj.DisplayRange = max(max(data_stat));
A=CGobj.ColumnLabels;
B=0;
close force 

plot(CGobj)

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir, subDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['Clustergram on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

% Create clusters

T2 = cluster(Z, 'cutoff', dendrocut, 'criterion', 'distance');

uniqueT2 = unique(T2); 

uniqueT2 = unique(T2); 
palette_couleurs = parula(numel(uniqueT2));

for i = 1:numel(T2)
    indice_cluster = T2(i);
    couleur_cluster = palette_couleurs(indice_cluster, :);
    T2(i, 2:4) = couleur_cluster; 
end

% PCA with clusters

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

scatter(zscores(:,1),zscores(:,2),50,T2(:,2:4),'filled','markeredgecolor',[0.5 0.5 0.5])
[p,q]=size(pc);

hold on

if numel(group) <= 15

[p,q]=size(pc);
hold on
for i=1:p
    if ((pc(i,1)^2+pc(i,2)^2)^0.5)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

elseif numel(group) > 15

[p,q]=size(pc);
hold on

distances = sqrt(pc(:,1).^2 + pc(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

end

xlabel(['PCA 1 ' num2str(fix(pcvars(1))) '% variance'],'fontsize',10);
ylabel(['PCA 2 ' num2str(fix(pcvars(2))) '% variance'],'fontsize',10);

% Clusters on the map 

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);
hold on;
for i = 1:m
    m_scatter(lon(i),lat(i), 50, T2(i,2:4), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
end

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCA and Map Clusters on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% PCoA plot with clusters

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

scatter(scores(:,1), scores(:,2), 50, T2(:,2:4), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
%text(scores(:,1),scores(:,2),sampleid,'HorizontalAlignment','center', 'VerticalAlignment','middle','Color','k','Interpreter','tex'); 
    
hold on 

if numel(group) <= 15

       n = size(vec,1); % get # of variables
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
   end

elseif numel(group) > 15

           n = size(vec,1); % get # of variables

distances = sqrt(vec(:,1).^2 + vec(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for j=topIndices'
    if distances(j)
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
    end
end

end

axis1 = sprintf('%2.2f',expl(1,1));
axis2 = sprintf('%2.2f',expl(2,1));
xlabel(['PCoA 1 (' num2str(axis1) ' % variance)'],'fontsize',10);
ylabel(['PCoA 2 (' num2str(axis2) ' % variance)'],'fontsize',10);

% Clusters on the map 

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);
hold on;
for i = 1:m
    m_scatter(lon(i),lat(i), 50, T2(i,2:4), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
end

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir, subDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCoA and Map Clusters on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% CA with clusters

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

subplot('Position', [0.05 0.1 0.3 0.8])

% distances
distances_X = sqrt(X(:,1).^2 + X(:,2).^2);
distances_Y = sqrt(Y(:,1).^2 + Y(:,2).^2);

if numel(group) <= 15

    %scatter(X(:,1), X(:,2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(:,1), X(:,2), group', 'Color', 'k', 'fontsize', 8);

    hold on
    for i=1:size(X,1)
    plot([0 X(i,1)],[0 X(i,2)],'k-')
    end

else

    [~, sortedIndices_X] = sort(distances_X, 'descend');
    topIndices_X = sortedIndices_X(1:round(0.3 * numel(distances_X)));
    
    %scatter(X(topIndices_X, 1), X(topIndices_X, 2), 30, 'filled', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    text(X(topIndices_X, 1), X(topIndices_X, 2), group(topIndices_X), 'Color', 'k', 'fontsize', 8);

    for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
    end

end

hold on;

scatter(Y(:,1), Y(:,2), 50, T2(:,2:4), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5]);

xlim([-max(abs([X(:,1); Y(:,1)])) * 1.1, max(abs([X(:,1); Y(:,1)])) * 1.1]);
ylim([-max(abs([X(:,2); Y(:,2)])) * 1.1, max(abs([X(:,2); Y(:,2)])) * 1.1]);

%text(Y(:,1) + .01, Y(:,2) + .1, sampleid', 'Color', 'k', 'fontsize', 8);

xlabel(['CA 1 (' num2str(fix(PI(1))) '% variance)']);
ylabel(['CA 2 (' num2str(fix(PI(2))) '% variance)']);

%Map with RGB colors

subplot('Position', [0.4 0.1 0.55 0.8])

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

[n,m] = size(data_stat);

[n,m] = size(data_stat);
hold on;
for i = 1:m
    m_scatter(lon(i),lat(i), 50, T2(i,2:4), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);
end

%Export and save? 
if doSave
    outputDir = fullfile(parentDir, subDir);
    mkdir(outputDir);
    
    filename = fullfile(outputDir, sprintf(['CA and Map Clusters on ' char(dataACP_name) ' ' list3{normalisation} '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

%% ── 12.6 NBSS spectra and slopes ────────────────────────────────────────────

if ~exist('ESD', 'var')
    ESD = base(1).tot.ESDvector;
end

% Truncated NBSS, all stations
figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');
for i = 1:numel(base)
    d = base(i).regroupped.NBSSliving_trunked;
    valid = d > 0 & isfinite(d);
    if any(valid)
        hold on
        plot(ESD(valid), d(valid), 'o', 'MarkerFaceColor', [0.9 0.9 0.9], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    end
end
set(gca, 'YScale', 'log', 'XScale', 'log');
xticks([0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000]);
xticklabels({'0.1', '0.2', '0.5', '1', '2', '5', '10', '20', '50', '100', '200', '500', '1000', '2000', '5000', '10000', '20000', '50000', '100000'});

xlabel('size on ESD (µm)', 'fontsize', 10);
ylabel({'NBSS (mm^3/mm^3/m^3)'; 'equivalent abundance'}, 'fontsize', 10);
title('NBSS of all stations', 'fontsize', 10);

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf('NBSS Spectra trunked.pdf'));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf

% NBSS slopes (reuse slope_values from EcoTaxa_Tool_Box.m section 11 if already computed)
if ~exist('slope_values', 'var')
    slope_values = arrayfun(@(b) b.regroupped.NBSSliving_slope, base);
end

% NBSS slope on PCA

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

scatter(zscores(:,1), zscores(:,2), 50, slope_values, 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5]); 

[p,q]=size(pc);

hold on

if numel(group) <= 15

[p,q]=size(pc);
hold on
for i=1:p
    if ((pc(i,1)^2+pc(i,2)^2)^0.5)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

elseif numel(group) > 15

[p,q]=size(pc);
hold on

distances = sqrt(pc(:,1).^2 + pc(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for i=topIndices'
    if distances(i)
    plot([0 pc(i,1)],[0 pc(i,2)],'k-')
    text(pc(i,1),pc(i,2),group(i),'fontsize',8);
    end
end

end

xlabel(['PCA 1 ' num2str(fix(pcvars(1))) '% variance'],'fontsize',10);
ylabel(['PCA 2 ' num2str(fix(pcvars(2))) '% variance'],'fontsize',10);

color = slandarerCM(6).Colors{1,9};  % spectral (change manually)
colormap(color);
colormap(flipud(color));
c2 = colorbar;

c=colorbar;
title(c,'NBSS Slope','fontweight','bold');

%annotation('ellipse', [0.65 0.025 0.01 0.015], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5, 'FaceColor', 'none');
%annotation('textbox', [0.66 0.02 0.6 0.04], 'String', '< 5 points on NBSS, slope not reliable', 'EdgeColor', 'none', 'FontSize', 8, 'Color', 'k');

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCA with NBSS slopes on ' char(dataACP_name) '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');

end

pause(1)
close gcf

% NBSS slope on PCoA

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

scatter(scores(:,1), scores(:,2), 50, slope_values, 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

hold on 

if numel(group) <= 10

       n = size(vec,1); % get # of variables
   for j = 1:n
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
   end

elseif numel(group) > 10

           n = size(vec,1); % get # of variables

distances = sqrt(vec(:,1).^2 + vec(:,2).^2);
[~, sortedIndices] = sort(distances, 'descend');
topIndices = sortedIndices(1:round(0.3 * numel(distances)));

for j=topIndices'
    if distances(j)
      % Plot vectors:
      f_arrow([0 0],[vec(j,1) vec(j,2)],'size',0.125*0.75,'angle',20,'Color','k')
      
      % Label vectors:
      h = text(vec(j,1)*delta(j),vec(j,2)*delta(j),group(j),'fontsize',8);
      set(h,'FontSize',8,'HorizontalAlignment','center','Color','k','Interpreter','tex');
    end
end

end

color = slandarerCM(6).Colors{1,9};  % spectral (change manually)
colormap(color);
colormap(flipud(color));
c2 = colorbar;

axis1 = sprintf('%2.2f',expl(1,1));
axis2 = sprintf('%2.2f',expl(2,1));
xlabel(['PCoA 1 (' num2str(axis1) ' % variance)'],'fontsize',10);
ylabel(['PCoA 2 (' num2str(axis2) ' % variance)'],'fontsize',10);

c=colorbar;
title(c,'NBSS Slope','fontweight','bold');

%annotation('ellipse', [0.65 0.025 0.01 0.015], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 0.5, 'FaceColor', 'none');
%annotation('textbox', [0.66 0.02 0.6 0.04], 'String', '< 5 points on NBSS, slope not reliable', 'EdgeColor', 'none', 'FontSize', 8, 'Color', 'k');

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir);
            mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['PCoA with NBSS slopes on ' char(dataACP_name) '.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end
pause(1)
close gcf

% NBSS slope on map

figure('numbertitle','off','Position',[10 50 scrsz(3) scrsz(4)-150],'Color','white');

m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
m_coast('patch',landcolor,'edgecolor',[0.6 0.6 0.6]);
hold on 
m_grid('fontsize',8,'linestyle','none','meridianlabel','off','parallellabel','off','xtick',[],'ytick',[]);

hold on;
m_scatter(lon(:), lat(:), 50, slope_values(:), 'filled', 'markeredgecolor', [0.5 0.5 0.5]);

color = slandarerCM(6).Colors{1,9};  % spectral (change manually)
colormap(color);
colormap(flipud(color));
c2 = colorbar;

c=colorbar;
title(c,'NBSS Slope','fontweight','bold');
%annotation('textbox', [0.7 0.05 0.2 0.1], 'String', 'empty if points nb on NBSS are < 5 : slope not reliable', 'EdgeColor', 'none', 'FontSize', 8, 'Color', 'k');

%Export and save? 
if doSave==false;
elseif doSave==true;

        outputDir = fullfile(parentDir);
        mkdir(outputDir);

    filename = fullfile(outputDir, sprintf(['Map of NBSS Slope.pdf']));
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
end

pause(1)
close gcf