function [decomp] = metisDecomp(fileName,numRegions,uval)
    consoleCMD = strcat(['gpmetis ','-ptype=kway ','-contig ','-objtype=cut ','-ncuts=1 ','-ctype=shem ',...
                    strcat(['-ufactor=',num2str(uval)]),' ',fileName,' ',num2str(numRegions)]);
    system(consoleCMD);
    resultName = strcat(fileName,'.part.',num2str(numRegions));
    
    fid = fopen(resultName,'r');
    decomp = fscanf(fid,'%i');
    fclose(fid);
    % system(strcat(['rm ',resultName]));
    system(strcat(['rm ',fileName]));
end