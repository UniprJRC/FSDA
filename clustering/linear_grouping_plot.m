function fh = linear_grouping_plot(X,hp,out)
% This function is called by lga and rlga and it is not intended to be called directly

% The peculiarity is the computation of hyperplanes for ortogonal
% regression. Arguments are:
% X = data
% hp = hyperplane
% out = output structure returned by function lga or rlga

% Copyright 2008-2021.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

k         = out.k;
alpha     = out.alpha;
plot_type = out.class;

[n,d]=size(X);

% this is just for rotating colors in the plots
clrdef = 'bkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcrbkmgcr';
symdef = '+*d^v><phos+*d^v><phos+*d^v><phos+*d^v><phos';

if d==2
    % initialize figure
    fh = figure('Name',[plot_type , ' plot'],'NumberTitle','off','Visible','on');
    gca(fh);
    hold on;
    
    % control of the axis limits
    xmin = min(X(:,1)); xmax = max(X(:,1));
    ymin = min(X(:,2)); ymax = max(X(:,2));
    deltax = (xmax - xmin) / 10;
    deltay = (ymax - ymin) / 10;
    
    xlim([xmin-deltax,xmax+deltax]);
    ylim([ymin-deltay,ymax+deltay]);
    
    % al contains the axes limits
    al = axis';
    
    % Plot the groups and the hyperplanes
    hold('on')
    for i=1:k
        group_label = ['Group ' num2str(i)];
        ucg = find(out.cluster==i);
        
        plot(X(ucg,1),X(ucg,2),'.w','DisplayName',[group_label ' (' num2str(length(ucg)) ' units)']);
        text(X(ucg,1),X(ucg,2),num2str(i*ones(length(ucg),1)),...
            'DisplayName',[group_label ' (' num2str(length(ucg)) ' units)'] , ...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'Color',clrdef(i));
        
        a =  hp(i, 3)/hp(i, 2);
        b = -hp(i, 1)/hp(i, 2);
        plot(al(1:2),a+b*al(1:2),'DisplayName',[group_label ' fit' ],'Color',clrdef(i));
    end
    
    % Plot the outliers (trimmed points), for rlga only
    if strcmp(plot_type , 'rlga')
        ucg = find(out.cluster==0);
        plot(X(ucg,1),X(ucg,2),'o','color','r','MarkerSize',8,...
            'DisplayName',['Trimmed units (' num2str(length(ucg)) ')']);
    end
    
    % Position the legends and make them clickable.
    lh=legend('show');
    legstr = get(lh,'String');
    clickableMultiLegend(legstr,'FontSize',14,'Location','northwest');
    
    axis('manual');
    
else
    % this is when p>2
    
    % axis labels
    nameY = cellstr([repmat('X',d,1) , num2str((1:d)')]);
    nameY = nameY';
    plo=struct;
    plo.nameY=nameY;
    plo.sym = [symdef(1:k) , 'o' ];
    plo.clr = [clrdef(1:k) , 'r' ];
    
    % group names in the legend
    group = cell(n,1);
    group(out.cluster==0) = {'Trimmed units'};
    for iii = 1:k
        group(out.cluster==iii) = {['Group ' num2str(iii)]};
    end
    
    spmplot(X,group,plo,'hist');
    
end

% add a title
if strcmp(plot_type , 'rlga')
    title([plot_type , ' plot' ' - ' '$\alpha = ' ...
        num2str(alpha) ' \; ROSS = ' num2str(out.ROSS) '$' ],...
        'interpreter' , 'LaTex', 'fontsize' , 18);
else
    title([plot_type , ' plot' ' - ' ' $\; ROSS = ' num2str(out.ROSS) '$' ],...
        'interpreter' , 'LaTex', 'fontsize' , 18);
end

end