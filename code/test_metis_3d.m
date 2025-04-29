% The following code assumes that you are in the same folder as this file
clear all;
close all;
restoredefaultpath();
addpath(genpath('../'));

maxNumCompThreads(15);

cm = readlines('cm.csv');
cm = cm(1:end-1);

degrees = 1:3;
subdivs = round(2.^([1:0.5:3]));
subdoms = 6:13;

%% Test parameters
runMETIS = false; % True requires a METIS installation
boolPlots = false;
view_a = -119.7638;
view_b = 26.3902;
plotSizeX = 24;
plotSizeY = 24;

%% Compute patch decomposition
sideLen = 3;
nrbarr = cnstrct_cube_nrbarr(sideLen*[1,1,1],[3,3,3]);

%% Problem Specification
% Physical parameters
mp_strct.nu  = @(x, y, z, ind) ones(size(x));

% Source and boundary terms
mp_strct.f = @(x, y, z, ind) cat(1,reshape(3.*cos(y).*cos(z).*sin(x),[1,size(x)]),...
                                   reshape(-6.*cos(x).*cos(z).*sin(y),[1,size(x)]),...
                                   reshape(3.*cos(x).*cos(y).*sin(z),[1,size(x)]));

% Exact solution (optional)
mp_strct.sol = @(x, y, z) cat(1,reshape(cos(y).*cos(z).*sin(x),[1,size(x)]),...
                                reshape(-2.*cos(x).*cos(z).*sin(y),[1,size(x)]),...
                                reshape(cos(x).*cos(y).*sin(z),[1,size(x)]));
mp_strct.sol_curl = @(x, y, z) cat(1,reshape(-3.*cos(x).*sin(y).*sin(z),[1,size(x)]),...
                                     zeros([1,size(x)]),...
                                     reshape(3.*cos(z).*sin(x).*sin(y),[1,size(x)]));

mp_strct.dir_bnd_func = @(x,y,z,iside) mp_strct.sol(x,y,z);

%% Export QoI
exportVals = numel(degrees)*numel(subdivs)*numel(subdoms);
err_curl = zeros(exportVals,1);
condEst = zeros(exportVals,1);
iter = zeros(exportVals,1);
numPrim = zeros(exportVals,1);
exportDeg = zeros(exportVals,1);
exportDivs = zeros(exportVals,1);
exportSubs = zeros(exportVals,1);
exportNumPatches = zeros(exportVals,1);


for dom=1:numel(subdoms)
    %% Compute decomposition with METIS
    numSub = subdoms(dom);
    if runMETIS
        graphName = strcat('cube',num2str(numel(nrbarr)),'.txt');
        mp2graphFile(nrbarr,graphName,[]);
        if numSub>2
            uval = 1;
        else
            uval = 50;
        end
        decomp = metisDecomp(graphName,numSub,uval);
    else
        fid = fopen(strcat('partitions/cube',num2str(numel(nrbarr)),'.txt.part.',num2str(numSub)),'r');
        decomp = fscanf(fid,'%i');
    end

    if ~all(ismember(0:numSub-1,decomp))
        error('METIS gave not the correct number of subdomains!')
    end

    %% Compute offset
    offset = {};
    for i=1:numSub
        subMean = [0,0,0]; subPatches = nrbarr(decomp+1==i);
        for j=1:numel(subPatches)
            subMean = subMean + nrbeval(subPatches(j),[0.5,0.5,0.5]);
        end
        subMean = subMean./numel(subPatches);
        offset{i} = 3*(0.5*[1.5,1.5,1.5] + 0.5*subMean);
    end
    
    %% Decomposition containing a subdomain with hole
    % decomp = [0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,2,2,2,2,2,2,2,2,2];

    %% Plot decomposition
    if boolPlots
        figure('Units','centimeters','Position',[0 0 plotSizeX plotSizeY]);
        clf()
        view(view_a,view_b);
        hold on;
        for i=1:numel(nrbarr)
            mynrbplot(nrbtform(nrbarr(i),vectrans(offset{decomp(i)+1})),[0,0,0],cm(decomp(i)+1));
        end
        hold off;
        axis off;
        % savefig(strcat('partitions/cube',num2str(numel(nrbarr)),'_part',num2str(numSub),'.fig'))
        % exportgraphics(gcf(),strcat(['decomp',num2str(numSub,'%02i'),'.pdf']),'ContentType','vector');
    end
    
    %% Get data for global problem
    [~,boundaries_glo,interfaces_glo,~,boundary_interfaces_glo] = mp_geo_load(nrbarr);
    
    %% Compute interface information
    subdomainCell = {};
    
    for iSub=0:numSub-1
        [geo,bounds,interfaces,~,b_intfcs] = mp_geo_load(nrbarr(decomp==iSub));
        
        mp_strct.patch_arr(iSub+1).geo = geo;
        mp_strct.patch_arr(iSub+1).bounds = bounds;
        mp_strct.patch_arr(iSub+1).intrfcs = interfaces;
        mp_strct.patch_arr(iSub+1).b_intrfcs = b_intfcs;
        
        subdomain.Geometry = geo;
        subdomain.Boundaries = bounds;
        subdomainCell{end+1} = subdomain;
    end
    [interfCell,boundarySidesCell] = getCouplingAndBoundaryInformation(subdomainCell,1e-6);
    
    for deg=1:numel(degrees)
        for sub=1:numel(subdivs)
    
            fprintf('\tStarting to simulate Nsub: %i, deg: %i, sub: %i\n',subdoms(dom),degrees(deg),subdivs(sub));
    
            %% Define canonical index
            canInd = 1 + (sub-1) + numel(subdivs)*(deg-1) + numel(subdivs)*numel(degrees)*(dom-1);
    
            %% Add discretization information
            degree = degrees(deg)*[1 1 1];
            nsub = subdivs(sub)*[1 1 1];
            nquad = degree+1;
            regularity = degree-1;
            [mp_strct] = add_basis_data(mp_strct,degree,nquad,regularity,nsub);
            
            %% Space and mesh construction of subdomains
            mshCell = {};
            spCell = {};
            spH1Cell = {};
            for iSub=1:numSub
                npatch = numel(mp_strct.patch_arr(iSub).geo);
                msh = cell (1, npatch);
                sp_hcurl  = cell (1, npatch);
                sp_h1  = cell (1, npatch);
                for iptc = 1:npatch
                    [knots, zeta] = ...
                        kntrefine (mp_strct.patch_arr(iSub).geo(iptc).nurbs.knots,...
                        mp_strct.patch_arr(iSub).nsub-1, mp_strct.patch_arr(iSub).degree,...
                        mp_strct.patch_arr(iSub).regularity);
                    [knots_hcurl, degree_hcurl] = knt_derham (knots, mp_strct.patch_arr(iSub).degree, 'Hcurl');
            
                    % Construct msh structure
                    rule      = msh_gauss_nodes (mp_strct.patch_arr(iSub).nquad);
                    [qn, qw]  = msh_set_quad_nodes (zeta, rule);
                    msh{iptc} = msh_cartesian (zeta, qn, qw, mp_strct.patch_arr(iSub).geo(iptc));
            
                    % Construct curl-space structure
                    scalar_spaces = cell (msh{iptc}.ndim, 1);
                    for idim = 1:msh{iptc}.ndim
                        scalar_spaces{idim} = sp_bspline (knots_hcurl{idim}, degree_hcurl{idim}, msh{iptc});
                    end
                    sp_hcurl{iptc} = sp_vector (scalar_spaces, msh{iptc}, 'curl-preserving');
                    clear scalar_spaces
            
                    % Construct grad-space structure
                    sp_h1{iptc} = sp_bspline (knots, mp_strct.patch_arr(iSub).degree, msh{iptc});
                end
                mshCell = cat(2,mshCell,msh);
                spCell = cat(2,spCell,sp_hcurl);
                spH1Cell = cat(2,spH1Cell,sp_h1);
            
                mp_strct.patch_arr(iSub).msh = msh_multipatch (msh, mp_strct.patch_arr(iSub).bounds);
                mp_strct.patch_arr(iSub).space = sp_multipatch (sp_hcurl, mp_strct.patch_arr(iSub).msh, mp_strct.patch_arr(iSub).intrfcs, mp_strct.patch_arr(iSub).b_intrfcs);
                mp_strct.patch_arr(iSub).space_h1 = sp_multipatch (sp_h1, mp_strct.patch_arr(iSub).msh, mp_strct.patch_arr(iSub).intrfcs, mp_strct.patch_arr(iSub).b_intrfcs);
                clear msh sp_h1 sp_hcurl
            end
            
            %% Space and mesh construction for global problem
            mp_strct.msh_mp = msh_multipatch (mshCell, boundaries_glo);
            mp_strct.space_mp = sp_multipatch (spCell, mp_strct.msh_mp, interfaces_glo, boundary_interfaces_glo);
            mp_strct.space_mul_mp = sp_multipatch (spH1Cell, mp_strct.msh_mp, interfaces_glo, boundary_interfaces_glo);
            
            %% Assemble local matrices and vectors
            mp_strct.cumu_dofs = 0;
            for iSub=1:numSub
                mp_strct.cumu_dofs = [mp_strct.cumu_dofs,mp_strct.cumu_dofs(end)+mp_strct.patch_arr(iSub).space.ndof];
                mp_strct.patch_arr(iSub).Ai = op_curlu_curlv_mp(mp_strct.patch_arr(iSub).space,...
                    mp_strct.patch_arr(iSub).space,...
                    mp_strct.patch_arr(iSub).msh, @(x,y,z) mp_strct.nu(x,y,z,iSub));
                mp_strct.patch_arr(iSub).fi = op_f_v_mp(mp_strct.patch_arr(iSub).space,...
                    mp_strct.patch_arr(iSub).msh, @(x,y,z) mp_strct.f(x,y,z,iSub))';
            
                mp_strct.patch_arr(iSub).Bi = double.empty(0,mp_strct.patch_arr(iSub).space.ndof);
            
                [mp_strct.patch_arr(iSub).gDi,...
                    mp_strct.patch_arr(iSub).dir_dofs] = sp_drchlt_l2_proj(mp_strct.patch_arr(iSub).space,...
                    mp_strct.patch_arr(iSub).msh, mp_strct.dir_bnd_func,...
                    boundarySidesCell{iSub});
            
                mp_strct.patch_arr(iSub).dep_dofs = [];
                mp_strct.patch_arr(iSub).idep_dofs = [];
            end
            
            %% Assemble coupling matrices
            for iSub = 1:numel(interfCell)
                intf = interfCell{iSub};
                dep_dofs = [];
                idep_dofs = [];
            
                for iBnd=1:numel(intf.sides1)
                    % Get information from both sides of mp-structure
                    i_patch = mp_strct.patch_arr(intf.reg1).bounds(intf.sides1(iBnd)).patches;
                    i_side = mp_strct.patch_arr(intf.reg1).bounds(intf.sides1(iBnd)).faces;
                    j_patch = mp_strct.patch_arr(intf.reg2).bounds(intf.sides2(iBnd)).patches;
                    j_side = mp_strct.patch_arr(intf.reg2).bounds(intf.sides2(iBnd)).faces;
                
                    % Introduce some checks if the interface-spaces are conforming
                    i_space = mp_strct.patch_arr(intf.reg1).space.sp_patch{i_patch}.boundary(i_side);
                    j_space = mp_strct.patch_arr(intf.reg2).space.sp_patch{j_patch}.boundary(j_side);
                
                    % save some index-stuff
                    i_gnum = mp_strct.patch_arr(intf.reg1).space.gnum{i_patch};
                    j_gnum = mp_strct.patch_arr(intf.reg2).space.gnum{j_patch};
                
                    % Save coupled DOFs
                    dep_dofs = union(dep_dofs,i_gnum(i_space.dofs),'stable');
                    idep_dofs = union(idep_dofs,j_gnum(j_space.dofs),'stable');
                end

                if numel(dep_dofs)~=numel(idep_dofs)
                    error('Something wrong with interface coupling!')
                end
            
                % Save dep und idep DOFs
                mp_strct.patch_arr(intf.reg1).dep_dofs = union(mp_strct.patch_arr(intf.reg1).dep_dofs,dep_dofs,'stable');
                mp_strct.patch_arr(intf.reg2).idep_dofs = union(mp_strct.patch_arr(intf.reg2).idep_dofs,idep_dofs,'stable');
            
                % Assemble coupling matrices
                B1 = spalloc(numel(dep_dofs),mp_strct.patch_arr(intf.reg1).space.ndof,0);
                B2 = spalloc(numel(idep_dofs),mp_strct.patch_arr(intf.reg2).space.ndof,0);
                
                B1(:,dep_dofs) = speye(numel(dep_dofs));
                B2(:,idep_dofs) = speye(numel(idep_dofs));
            
                mp_strct.patch_arr(intf.reg1).Bi = [mp_strct.patch_arr(intf.reg1).Bi;B1];
                mp_strct.patch_arr(intf.reg2).Bi = [mp_strct.patch_arr(intf.reg2).Bi;-B2];
                for iPatch = setdiff(1:numSub,[intf.reg1,intf.reg2])
                    mp_strct.patch_arr(iPatch).Bi = [mp_strct.patch_arr(iPatch).Bi;...
                        zeros([size(B1,1),mp_strct.patch_arr(iPatch).space.ndof])];
                end
            end
            
            %% Construct graphs on patch, subdomain and global level
            [mp_strct] = graphsWithPrio_METIS(mp_strct,[interfCell{:}],boundarySidesCell,decomp,false);

            %% Plot wire basket
            if boolPlots
                figure('Units','centimeters','Position',[0 0 plotSizeX plotSizeY]);
                clf()
                view(view_a,view_b);
                hold on;
                for i=1:numSub
                    marked = {[mp_strct.patch_arr(i).graph.Edges.IDs(mp_strct.patch_arr(i).graph.Edges.Weight==1)',...
                               mp_strct.patch_arr(i).graph.Edges.IDs(mp_strct.patch_arr(i).graph.Edges.Weight==2)'],...
                               mp_strct.patch_arr(i).graph.Edges.IDs(mp_strct.patch_arr(i).graph.Edges.Weight==3)'};
                    plot_graph_and_marked(mp_strct.patch_arr(i).graph,marked,[],offset{i},[],[],4,10,["#005aa9", "#9ac103"]);
                end
                hold off;
                axis off;
                exportgraphics(gcf(),strcat(['decomp',num2str(numSub,'%02i'),'_wirebasket.pdf']),'ContentType','vector');
            end

            %% Construct global tree and project on patch and subdomain level
            % Construct global tree
            T = minspantree(mp_strct.gloGraph,'Method','sparse');
            mp_strct.gloTree = T.Edges.IDs;
            
            wirebasketDOFs = mp_strct.gloGraph.Edges.IDs(mp_strct.gloGraph.Edges.Weight==2);
            mp_strct.gloPri = setdiff(wirebasketDOFs,mp_strct.gloTree);
            
            % Project special DOFs onto patch graphs
            patchTreeComps = cell(1,mp_strct.space_mp.npatch);
            patchPriComps = cell(1,mp_strct.space_mp.npatch);
            for i=1:mp_strct.space_mp.npatch
                gDOFs = mp_strct.space_mp.gnum{i};
                lDOFs = 1:numel(gDOFs);
                patchTreeComps{i} = lDOFs(ismember(gDOFs,mp_strct.gloTree));
                patchPriComps{i} = lDOFs(ismember(gDOFs,mp_strct.gloPri));
    
                for j=1:numel(lDOFs)
                    patchBool = mp_strct.patchGraphs{i}.Edges.IDs==j;
                    gloBool = mp_strct.gloGraph.Edges.IDs==gDOFs(j);
                    mp_strct.patchGraphs{i}.Edges.Weight(patchBool) = mp_strct.gloGraph.Edges.Weight(gloBool);
                end
            end
            
            % Project tree from patch level onto subdomain level
            for i=1:numSub
                patchIDs = find(decomp==i-1);
                
                treeComps = patchTreeComps(patchIDs);
                priComps = patchPriComps(patchIDs);
                subTreeComps = [];
                subPriComps = [];
                for j=1:numel(treeComps)
                    gDOFs = mp_strct.patch_arr(i).space.gnum{j};
                    subTreeComps = union(subTreeComps,gDOFs(treeComps{j}));
                    subPriComps = union(subPriComps,gDOFs(priComps{j}));
                    
     
                    for k=1:numel(gDOFs)
                        patchBool = mp_strct.patchGraphs{i}.Edges.IDs==k;
                        subBool = mp_strct.patch_arr(i).graph.Edges.IDs==gDOFs(k);
                        mp_strct.patch_arr(i).graph.Edges.Weight(subBool) = mp_strct.patchGraphs{patchIDs(j)}.Edges.Weight(patchBool);
                    end
                end
                mp_strct.patch_arr(i).tree = subTreeComps;
                mp_strct.patch_arr(i).pri_dofs = subPriComps;
    
    
            end
            
            %% Compute Cp
            BpCell = cellfun(@(B,pri) B(:,pri), {mp_strct.patch_arr.Bi}, {mp_strct.patch_arr.pri_dofs}, 'UniformOutput',false);
            Bp = [BpCell{:}];
            
            redBp = Bp(sum(abs(Bp),2)~=0,:);
            Cp = null(full(redBp));
            
            low = 1;
            for i=1:numSub
                upp = low - 1 + numel(mp_strct.patch_arr(i).pri_dofs);
                mp_strct.patch_arr(i).Cp = Cp(low:upp,:);
                low = upp + 1;
            end
    
    
            if boolPlots
                figure('Units','centimeters','Position',[0 0 plotSizeX plotSizeY]);
                clf()
                view(view_a,view_b);
                % title('Exploded View of Subdomains')
                hold on;
                for i=1:numSub
                    marked = {mp_strct.patch_arr(i).graph.Edges.IDs(mp_strct.patch_arr(i).graph.Edges.Weight==1),...
                        mp_strct.patch_arr(i).graph.Edges.IDs(mp_strct.patch_arr(i).graph.Edges.Weight==2),...
                        mp_strct.patch_arr(i).graph.Edges.IDs(mp_strct.patch_arr(i).graph.Edges.Weight==3),...
                        mp_strct.patch_arr(i).graph.Edges.IDs(mp_strct.patch_arr(i).graph.Edges.Weight==4),...
                        mp_strct.patch_arr(i).graph.Edges.IDs(mp_strct.patch_arr(i).graph.Edges.Weight==5)};
                    plot_graph_and_marked(mp_strct.patch_arr(i).graph,marked,[],offset{i},[],[],2,10,["#005aa9", "#ec6400", "#9ac103", "#710c84", "#000000"]);
                end
                hold off;
                axis off;
                exportgraphics(gcf(),strcat(['decomp',num2str(numSub,'%02i'),'_weights.pdf']),'ContentType','vector');
            end

            if boolPlots
                figure('Units','centimeters','Position',[0 0 plotSizeX plotSizeY]);
                clf()
                view(view_a,view_b);
                hold on;
                for i=1:numSub
                    marked = {mp_strct.patch_arr(i).tree, mp_strct.patch_arr(i).pri_dofs};
                    plot_graph_and_marked(mp_strct.patch_arr(i).graph,marked,[],offset{i},[],[],2,10,["#009D81", "#A60084"]);
                end
                hold off;
                axis off;
                exportgraphics(gcf(),strcat(['decomp',num2str(numSub,'%02i'),'_tree.pdf']),'ContentType','vector');
            end
            
            %% Solve system with ieti-dp
            [sol,multi,condEst(canInd),iter(canInd)] = metis_solver(mp_strct,1e-6);
            numPrim(canInd) = size(Cp,2);
            exportDeg(canInd) = degrees(deg);
            exportDivs(canInd) = subdivs(sub);
            exportNumPatches(canInd) = numel(decomp);
            exportSubs(canInd) = numSub;
            [~,~,err_curl(canInd)] = hcurl_error_computation(mp_strct,sol);
            fprintf('\n');
        end
    end
end

%% Export data
tab = table();
tab.deg = exportDeg;
tab.divs = exportDivs;
tab.patchs = exportNumPatches;
tab.subs = exportSubs;
tab.pri = numPrim;
tab.cond = condEst;
tab.iter = iter;
tab.err = err_curl;

writetable(tab,'data/generated_data.csv');