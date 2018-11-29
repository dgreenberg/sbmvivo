function fn = rfilename(oedb, fileinfo, algind, dataset, neuron, segment)
fn = [fileinfo.tmpdir 'r' num2str(algind) '_' oedb.segcode{dataset}{neuron}{segment} '.mat'];