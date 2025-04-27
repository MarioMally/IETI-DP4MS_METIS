function [subGraph,newLocGraphs] = loc2sub_graph(graphCell,mp_space,mp_mul_space,intfcs,bounds,conSub)
    
    gloNodes = [1:mp_mul_space.ndof]';
%     gloNodePrio = zeros(mp_mul_space.ndof,1);
    gloEdges = zeros(mp_space.ndof,1);
    gloEndNodes = zeros(mp_space.ndof,2);
    gloEdgePrio = zeros(mp_space.ndof,1);
    
    xCoor = zeros(numel(gloNodes),1);
    yCoor = zeros(numel(gloNodes),1);
    zCoor = zeros(numel(gloNodes),1);
%     geoInfo = [{},{}];

    newLocGraphs = cell(1,numel(graphCell));
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

        boolPatches = [intfcs.patches] == iPatch;
        intrfcFaces = [intfcs(boolPatches).faces];
        patchConSub = conSub(boolPatches);

        boundFaces = [bounds([bounds.patches]==iPatch).faces];
        locWeights = zeros(numel(ids),1);

        for jEdge=1:numel(ids)
            if isempty(G.Edges.geoInfo{jEdge})
                
                locWeights(jEdge) = 5; % Part of Interior
                gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);

            elseif numel(G.Edges.geoInfo{jEdge})==1

                geoInfo = G.Edges.geoInfo{jEdge};
                if ismember(geoInfo,boundFaces) || ismember(geoInfo,intrfcFaces)
                    locWeights(jEdge) = 4; % Part of facet
                    gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);
                else
                    locWeights(jEdge) = 5; % Part of Interior
                    gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);
                end
            elseif numel(G.Edges.geoInfo{jEdge})==2

                geoInfo = G.Edges.geoInfo{jEdge};
                if ismember(geoInfo(1),intrfcFaces) && ismember(geoInfo(2),boundFaces)
                    locWeights(jEdge) = 1; % Intersection Interface and Boundary
                    gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);
                elseif ismember(geoInfo(1),boundFaces) && ismember(geoInfo(2),intrfcFaces)
                    locWeights(jEdge) = 1; % Intersection Interface and Boundary
                    gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);
                elseif ismember(geoInfo(1),intrfcFaces) && ismember(geoInfo(2),intrfcFaces)
                    [~,ind1] = ismember(geoInfo(1),intrfcFaces);
                    [~,ind2] = ismember(geoInfo(2),intrfcFaces);
                    if patchConSub(ind1)==patchConSub(ind2) % check if the connected interface sides belong to the same subdomain
                        locWeights(jEdge) = 4; % Part of standard interface
                    else
                        locWeights(jEdge) = 2; % Cross-Edges
                    end
                    gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);
                elseif ismember(geoInfo(1),boundFaces) && ismember(geoInfo(2),boundFaces)
                    locWeights(jEdge) = 3; % Wire basket, Dirichlet edges
                    gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);
                elseif any(ismember(geoInfo,boundFaces)) || any(ismember(geoInfo,intrfcFaces))
                    locWeights(jEdge) = 4; % Part of facet
                    gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);
                else
                    locWeights(jEdge) = 5; % Part of interior
                    gloEdgePrio(gEdges(ids(jEdge))) = locWeights(jEdge);
                end

            else
                error('Something went wrong with geoInfo of Edges');
            end

            newNodeTable = table();
            newEdgeTable = table();
            
            newNodeTable.Name = G.Nodes.Name;
            newNodeTable.IDs = G.Nodes.IDs;
            newNodeTable.xCoor = G.Nodes.xCoor;
            newNodeTable.yCoor = G.Nodes.yCoor;
            newNodeTable.zCoor = G.Nodes.zCoor;

            newEdgeTable.EndNodes = G.Edges.EndNodes;
            newEdgeTable.IDs = G.Edges.IDs;
            newEdgeTable.Weight = locWeights;
            newLocGraphs{iPatch} = graph(newEdgeTable,newNodeTable);
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
    EdgeTable.Weight = gloEdgePrio;

    subGraph = graph(EdgeTable,NodeTable);

end