function mp2graphFile(nrbarr,fileName,matTypeLin)
    [intf,~] = nrbmultipatch(nrbarr);
    nodes1 = [intf.patch1];
    nodes2 = [intf.patch2];

    fid = fopen(fileName,'w');
    if ~isempty(matTypeLin)
        fprintf(fid,'%i %i 001\n',numel(nrbarr),numel(intf));
    else
        fprintf(fid,'%i %i\n',numel(nrbarr),numel(intf));
    end

    for i=1:numel(nrbarr)
        if ~isempty(matTypeLin)
            iIsLin = ismember(i,matTypeLin);
        end
        edges = union(nodes2(nodes1==i),nodes1(nodes2==i));
        for j=1:numel(edges)
            if ~isempty(matTypeLin)
                jIsLin = ismember(j,matTypeLin);
                if iIsLin && jIsLin
                    weight = 10000;
                elseif iIsLin || jIsLin
                    weight = 1;
                else
                    weight = 100;
                end
                fprintf(fid,'%i %i ',edges(j),weight);
            else
                fprintf(fid,'%i ',edges(j));
            end
        end
        fprintf(fid,'\n');
    end

    fclose(fid);

end