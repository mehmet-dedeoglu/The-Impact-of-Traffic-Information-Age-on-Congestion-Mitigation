function [G,G_ind,G_weight,G_cap,G_adj,M] = createGraph(w,choice,M)

% if nargin == 0
if choice == 1 
    G_ind = [2,3,7,49,49;1,3,10,49,49;1,2,4,49,49;3,5,49,49,49;4,6,48,49,49;...
        5,7,49,49,49;1,6,8,21,49;7,9,11,18,49;8,10,49,49,49;2,9,20,49,49;...
        8,12,13,18,49;11,14,49,49,49;11,14,49,49,49;12,13,15,16,49;...
        14,49,49,49,49;14,17,27,49,49;16,18,21,49,49;8,11,17,19,20;...
        18,20,49,49,49;10,18,19,41,49;7,17,22,49,49;21,23,49,49,49;...
        22,24,28,48,49;23,25,49,49,49;24,29,40,49,49;27,28,49,49,49;...
        16,26,29,49,49;23,26,36,49,49;25,27,30,49,49;29,31,32,35,38;...
        30,32,49,49,49;30,31,33,34,35;32,49,49,49,49;32,49,49,49,49;...
        30,32,36,37,49;28,35,49,49,49;35,38,39,43,49;30,37,42,44,49;...
        37,40,49,49,49;25,39,41,46,49;20,40,42,43,49;38,41,43,49,49;...
        37,41,42,49,49;38,47,49,49,49;46,49,49,49,49;40,45,47,49,49;...
        44,46,48,49,49;5,23,47,49,49];
    G_weight = ones(48,5);
    G_weighta = G_weight';
    G_weighta(G_ind'~=49)=w;
    G_weight = G_weighta';
%     
%     G_weight(1,2) = 3;
%     
    G_index = zeros(size(G_ind));
    for k = 1:48
        G_index(k,:) = (k-1)*49*ones(1,5)+G_ind(k,:);
    end
    G_adj = zeros(49);
    G_adj(G_index) = 1;
    G_adj = G_adj(1:48,1:48);
    M = sum(sum(G_adj));
    
    G_cap = zeros(49);
    G_cap(G_index) = 1000;
    G_cap = G_cap(1:48,1:48);
%     
%     weight(2) = 3;
%     weight(7) = 3;
%     figure;
    G_adj(G_adj==1) = G_adj(G_adj==1).*w';
    G = digraph(G_adj);
%     G_simple = graph(G_adj);
%     fig1 = plot(G_simple);
%     labeledge(fig1,1:numedges(G),1:numedges(G))


% figure(1)
% h = plot(G_simple,'EdgeLabel',G_simple.Edges.Weight);
% nl = h.NodeLabel;
% h.NodeLabel = '';
% h.EdgeLabel = '';
% xd = get(h, 'XData');
% yd = get(h, 'YData');
% text(xd+0.05, yd, nl, 'FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')


%     p = shortestpath(G,1,45);
%     highlight(fig1,'59','EdgeColor','r');
    
%     TR = shortestpathtree(G,1);
%     figure;plot(TR);
% end
elseif choice == 2
    G_ind = [2,7,9,9,9;1,3,8,9,9;2,4,9,9,9;3,5,9,9,9;4,6,8,9,9;...
        5,7,9,9,9;1,6,8,9,9;2,5,7,9,9];
    G_weight = ones(8,5);
    G_weighta = G_weight';
    G_weighta(G_ind'~=9)=w;
    G_weight = G_weighta';

    G_index = zeros(size(G_ind));
    for k = 1:8
        G_index(k,:) = (k-1)*9*ones(1,5)+G_ind(k,:);
    end
    G_adj = zeros(9);
    G_adj(G_index) = 1;
    G_adj = G_adj(1:8,1:8);
    M = sum(sum(G_adj));
    
    G_cap = zeros(49);
    G_cap(G_index) = 1000;
    G_cap = G_cap(1:48,1:48);

    G_adj(G_adj==1) = G_adj(G_adj==1).*w';
    G = digraph(G_adj);
    
elseif choice == 3
    G_ind = [2,3,24,24,24,24;1,4,5,6,7,8;1,9,10,12,27,24;2,13,24,24,24,24;...
        2,9,13,14,15,16;2,14,16,17,24,24;2,12,18,24,24,24;2,14,24,24,24,24;...
        3,5,19,24,24,24;3,20,24,24,24,24;3,20,21,24,24,24;3,7,14,19,22,24;...
        4,5,19,24,24,24;5,6,8,12,23,24;5,16,24,24,24,24;5,6,15,24,24,24;...
        6,19,24,24,24,24;7,23,24,24,24,24;9,12,13,17,22,24;10,11,24,24,24,24;...
        11,23,24,24,24,24;12,19,24,24,24,24;14,18,21,24,24,24];
    G_weight = ones(23,6);
    G_weighta = G_weight';
    G_weighta(G_ind'~=24)=w;
    G_weight = G_weighta';

    G_index = zeros(size(G_ind));
    for k = 1:23
        G_index(k,:) = (k-1)*24*ones(1,6)+G_ind(k,:);
    end
    G_adj = zeros(24);
    G_adj(G_index) = 1;
    G_adj = G_adj(1:23,1:23);
    M = sum(sum(G_adj));
    
    G_cap = zeros(24);
    G_cap(G_index) = 1000;
    G_cap = G_cap(1:23,1:23);

    G_adj(G_adj==1) = G_adj(G_adj==1).*w';
    G = digraph(G_adj);
    
elseif choice == 4
    G_ind = [3,4,6,7,9,10,14,16,16,16;14,16,16,16,16,16,16,16,16,16;1,4,7,9,10,11,13,14,16,16;...
        1,3,7,9,10,14,16,16,16,16;6,16,16,16,16,16,16,16,16,16;1,5,9,11,13,16,16,16,16,16;...
        1,3,4,9,10,13,14,16,16,16;10,16,16,16,16,16,16,16,16,16;1,3,4,6,7,10,12,13,14,15;...
        1,3,4,7,8,9,13,14,16,16;3,6,16,16,16,16,16,16,16,16;9,16,16,16,16,16,16,16,16,16;...
        3,6,7,9,10,14,16,16,16,16;1,2,3,4,7,9,10,13,16,16;9,16,16,16,16,16,16,16,16,16];
    G_weight = ones(15,10);
    G_weighta = G_weight';
    G_weighta(G_ind'~=16)=w;
    G_weight = G_weighta';

    G_index = zeros(size(G_ind));
    for k = 1:15
        G_index(k,:) = (k-1)*16*ones(1,10)+G_ind(k,:);
    end
    G_adj = zeros(16);
    G_adj(G_index) = 1;
    G_adj = G_adj(1:15,1:15);
    M = sum(sum(G_adj));
    
    G_cap = zeros(16);
    G_cap(G_index) = 1000;
    G_cap = G_cap(1:15,1:15);

    G_adj(G_adj==1) = G_adj(G_adj==1).*w';
    G = digraph(G_adj);

elseif choice == 5%Not complete
    G_ind = [2,3,4,5,6,7,8,9,10,11,15;1,15,15,15,15,15,15,15,15,15,15;...
        1,15,15,15,15,15,15,15,15,15,15;1,15,15,15,15,15,15,15,15,15,15;...
        1,15,15,15,15,15,15,15,15,15,15;1,7,15,15,15,15,15,15,15,15,15;...
        1,6,15,15,15,15,15,15,15,15,15;1,9,10,12,15,15,15,15,15,15,15;...
        1,8,10,11,12,15,15,15,15,15,15;1,8,9,11,12,13,14,15,15,15,15;...
        1,9,10,13,15,15,15,15,15,15,15;8,9,10,13,15,15,15,15,15,15,15;...
        10,11,12,15,15,15,15,15,15,15,15;10,15,15,15,15,15,15,15,15,15,15];
    G_weight = ones(14,10);
    G_weighta = G_weight';
    G_weighta(G_ind'~=15)=w;
    G_weight = G_weighta';

    G_index = zeros(size(G_ind));
    for k = 1:14
        G_index(k,:) = (k-1)*15*ones(1,10)+G_ind(k,:);
    end
    G_adj = zeros(15);
    G_adj(G_index) = 1;
    G_adj = G_adj(1:14,1:14);
    M = sum(sum(G_adj));
    
    G_cap = zeros(15);
    G_cap(G_index) = 1000;
    G_cap = G_cap(1:14,1:14);

    G_adj(G_adj==1) = G_adj(G_adj==1).*w';
    G = digraph(G_adj);

end

end