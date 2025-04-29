clear all

paper_data = readmatrix("paper_data.csv","Delimiter",",");
gener_data = readmatrix("generated_data.csv","Delimiter",",");

assert(all(size(paper_data) == size(gener_data)));

for i=1:size(paper_data,2)
    rel_err = norm(paper_data(:,i)-gener_data(:,i))./norm(paper_data(:,i));
    fprintf('Column %i with rel. error %0.3d\n',i,rel_err);
    assert(rel_err<=1e-10);
end