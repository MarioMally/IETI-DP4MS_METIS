function [gloGraph] = loc2glo_graph(graphCell,mp_space,mp_mul_space,offset)
    
    gloNodes = [1:mp_mul_space.ndof]';
%     gloNodePrio = zeros(mp_mul_space.ndof,1);
    gloEdges = zeros(mp_space.ndof,1);
    gloEndNodes = zeros(mp_space.ndof,2);
    gloWeight = inf*ones(mp_space.ndof,1);
    
    xCoor = inf*ones(numel(gloNodes),1);
    yCoor = inf*ones(numel(gloNodes),1);
    zCoor = inf*ones(numel(gloNodes),1);
%     geoInfo = [{},{}];

    for iPatch=1:numel(graphCell)
        G = graphCell{iPatch};
        % Node stuff
        gNodes = mp_mul_space.gnum{iPatch};
        xCoor(gNodes) = G.Nodes.xCoor;
        yCoor(gNodes) = G.Nodes.yCoor;
        zCoor(gNodes) = G.Nodes.zCoor;
        

        % Edge stuff
        gEdges = mp_space.gnum{iPatch};
        ids = G.Edges.IDs;
        endNodes = str2double(G.Edges.EndNodes);

        gloEndNodes(gEdges(ids),:) = gNodes(endNodes);
        gloEdges(gEdges(ids)) =  gEdges(ids);

        for j=1:numel(ids)
            id = gEdges(ids(j));
            % Check for parabolic problem to identify intersection between conducting and nonconducting regions
            if (gloWeight(id)==2 &&  G.Edges.Weight(j)==offset + 2) || (gloWeight(id)==offset+2 &&  G.Edges.Weight(j)==2)
                gloWeight(id) = 1;
            elseif (gloWeight(id)==4 &&  G.Edges.Weight(j)==offset + 4) || (gloWeight(id)==offset + 4 &&  G.Edges.Weight(j)==4)
                gloWeight(id) = 3;
            elseif (gloWeight(id)==7 &&  G.Edges.Weight(j)==offset + 7) || (gloWeight(id)==offset + 7 &&  G.Edges.Weight(j)==7)
                gloWeight(id) = 6;
            else
                gloWeight(id) = min([gloWeight(id),G.Edges.Weight(j)]);
            end
        end
    end

    % boolW1 = gloWeight==1;
    % boolW2 = gloWeight==2;
    % boolW3 = gloWeight==3;
    % boolW4 = gloWeight==4;
    % boolW6 = gloWeight==6;
    % boolW7 = gloWeight==7;
    % gloWeight(boolW1) = 2;
    % gloWeight(boolW2) = 1;
    % gloWeight(boolW3) = 4;
    % gloWeight(boolW4) = 3;
    % gloWeight(boolW6) = 7;
    % gloWeight(boolW7) = 6;
    
    %% Set up graph from tables
    NodeTable = table();
    EdgeTable = table();

    NodeTable.Name = cellfun(@(x) int2str(x), num2cell(gloNodes),'UniformOutput',false);
    NodeTable.IDs = gloNodes;
    NodeTable.xCoor = xCoor;
    NodeTable.yCoor = yCoor;
    NodeTable.zCoor = zCoor;
%     NodeTable.geoInfo = geoInfo';
%     NodeTable.Weight = gloNodePrio;

    EdgeTable.EndNodes = cellfun(@(x) int2str(x), num2cell(gloEndNodes),'UniformOutput',false);
    EdgeTable.IDs = gloEdges;
    EdgeTable.Weight = gloWeight;

    gloGraph = graph(EdgeTable,NodeTable);

end