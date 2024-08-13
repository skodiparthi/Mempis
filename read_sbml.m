% model = readCbModel;

IndexC = strfind(model.genes, 'G5_RASTassembly.CDS.5357');
Index = find(not(cellfun('isempty',IndexC)))