function [mp_strct] = graphsWithPrio_METIS(mp_strct,intrfc,boundarySidesCell,decomp,plotting)
    
    if plotting

        figure(1)
        clf();
        
        nrbcell = cell(1,numel(decomp));

        ax2 = subplot(2,2,2);
        subtitle('Subdomain Graphs');
    end

    patchGraphs = cell(1,numel(decomp));
    for iSub=1:numel(mp_strct.patch_arr)

        % Relevant space and offset
        space_h1 = mp_strct.patch_arr(iSub).space_h1;
        space = mp_strct.patch_arr(iSub).space;
        bounds = mp_strct.patch_arr(iSub).bounds(boundarySidesCell{iSub});

        sides = [];
        conSub = [];
        for i=1:numel(intrfc)
            sub1 = intrfc(i).reg1;
            sub2 = intrfc(i).reg2;
            if iSub == sub1
                sides = [sides, intrfc(i).sides1];
                conSub = [conSub, sub2*ones(size(intrfc(i).sides1))];
            elseif iSub == sub2
                sides = [sides, intrfc(i).sides2];
                conSub = [conSub, sub1*ones(size(intrfc(i).sides2))];
            end
        end
        intrfcs = mp_strct.patch_arr(iSub).bounds(sides);
        
        %% Patch graphs
        subPatchGraphs = cell(1,space_h1.npatch);
        ids = find(ismember(decomp,iSub-1));
        for j=1:space_h1.npatch
            % Some parameter to support plotting of multiple graphs in one plot
            nrb = mp_strct.patch_arr(iSub).geo(j).nurbs;
            if plotting
                nrbcell{ids(j)} = nrb;
            end
            % patch graph construction
            subPatchGraphs{j} = get_system_graph2(space_h1.sp_patch{j},0,0,nrb);
        end

        %% Put together patch graphs to subdomain graph
        [subGraph,subPatchGraphs] = loc2sub_graph(subPatchGraphs,space,space_h1,intrfcs,bounds,conSub);
        patchGraphs(ids) = subPatchGraphs;
        clear subPatchGraphs

        %% Some plotting
        if plotting
            e1 = find(subGraph.Edges.Weight==1 | subGraph.Edges.Weight==2);
            e2 = find(subGraph.Edges.Weight==3);
            e3 = find(subGraph.Edges.Weight==4);
            markedSub = {subGraph.Edges.IDs(e1),subGraph.Edges.IDs(e2)};%,subGraph.Edges.IDs(e3)};
            
            hold on;
            plot_graph_and_marked(subGraph,markedSub,[],[0,0,0],[],[],2,5);
        end
        
        mp_strct.patch_arr(iSub).graph = subGraph;
    end


    %% Plot patch graphs with prio
    if plotting
        ax3 = subplot(2,2,3);
        subtitle('Patch Graphs');
        for j=1:numel(patchGraphs)
            e1 = find(patchGraphs{j}.Edges.Weight==1 | patchGraphs{j}.Edges.Weight==2);
            e2 = find(patchGraphs{j}.Edges.Weight==3);
            e3 = find(patchGraphs{j}.Edges.Weight==4);
            markedLoc = {patchGraphs{j}.Edges.IDs(e1),patchGraphs{j}.Edges.IDs(e2)};%,patchGraphs{j}.Edges.IDs(e3)};
    
            hold on;
            plot_graph_and_marked(patchGraphs{j},markedLoc,[],[0,0,0],[],[],2,5);
        end
    end
    
    %% Compute global graph
    gloGraph = loc2glo_graph(patchGraphs,mp_strct.space_mp,mp_strct.space_mul_mp,0);

    if plotting
        ddcolormap = ["#34a8dc", "#a25264", "#85c79c", "#cda19a", "#abd533", "#c85221", "#8886e9", "#12d388", "#4b7359", "#eb67f9", "#f5bd5e", "#cc156d", "#ab7b05", "#63e118", "#fd2c3b"];
        ax1 = subplot(2,2,1);
        subtitle('Decomposition');
        hold on;
        for i=1:numel(nrbcell)
            mynrbplot(nrbcell{i},[0,0,0],ddcolormap(decomp(i)+1));
        end
        hold off;
        axis off;

        ax4 = subplot(2,2,4);
        subtitle('Global Graph');
    
        e1 = find(gloGraph.Edges.Weight==1 | gloGraph.Edges.Weight==2);
        e2 = find(gloGraph.Edges.Weight==3);
        e3 = find(gloGraph.Edges.Weight==4);
        markedGlo = {gloGraph.Edges.IDs(e1),gloGraph.Edges.IDs(e2)};%,gloGraph.Edges.IDs(e3)};
        hold on;
        plot_graph_and_marked(gloGraph,markedGlo,[],[0,0,0],[],[],2,5);

        Link = linkprop([ax1, ax2, ax3, ax4],{'CameraUpVector', 'Projection', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);
        ax1.PlotBoxAspectRatioMode = 'auto';
        ax1.DataAspectRatio = [1,1,1];
    end

    mp_strct.gloGraph = gloGraph;
    mp_strct.patchGraphs = patchGraphs;
end