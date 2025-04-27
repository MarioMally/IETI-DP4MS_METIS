function [interfCell,boundarySidesCell] = getCouplingAndBoundaryInformation(regionCell,tol)
    
    interfCell = {};
    boundarySidesCell = {};
    for i=1:numel(regionCell)
        boundarySidesCell{i} = 1:numel(regionCell{i}.Boundaries);
    end

    for i=1:numel(regionCell)-1
        for j=i+1:numel(regionCell)
%             fprintf('%i, %i, ',i,j);
            iRegion = regionCell{i};
            jRegion = regionCell{j};
            iBoundaries = iRegion.Boundaries;
            jBoundaries = jRegion.Boundaries;
    
            isides = [];
            jsides = [];
            ornt1 = [];
            ornt2 = [];
            
            % Naively Check Every boundary (more clever maybe later)
            for k=1:numel(iRegion.Geometry)
                for l=1:numel(jRegion.Geometry)

                    [intf,~] = nrbmultipatch([iRegion.Geometry(k).nurbs,jRegion.Geometry(l).nurbs],tol);
                    
                    if ~isempty(intf)
                        isides(end+1) = find(([iBoundaries.patches]==k) & ([iBoundaries.faces]==intf.side1));
                        jsides(end+1) = find(([jBoundaries.patches]==l) & ([jBoundaries.faces]==intf.side2));
                        ornt1(end+1) = intf.ornt1;
                        ornt2(end+1) = intf.ornt2;
                    end
                end
            end
    
            if ~isempty(isides) && ~isempty(jsides)
%                 fprintf('true\n')
                % Save interface data
                interf.reg1 = i;
                interf.reg2 = j;
                interf.sides1 = isides;
                interf.sides2 = jsides;
                interf.ornt1 = ornt1;
                interf.ornt2 = ornt2;
                interfCell{end+1} = interf;
    
                %remove boundary data
                boundarySidesCell{i} = setdiff(boundarySidesCell{i},isides);
                boundarySidesCell{j} = setdiff(boundarySidesCell{j},jsides);
            else
%                 fprintf('false\n')
            end
        end
    end
end