function [gloGraph] = loc2glo_graph(graphCell,mp_space,mp_mul_space)
    
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
            gloWeight(id) = min([gloWeight(id),G.Edges.Weight(j)]);
        end

    end
    
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