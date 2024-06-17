import stereo as st
import sys
data_path=sys.argv[1]
data = st.io.read_gem(file_path=data_path, bin_type='bins', bin_size=40, is_sparse=True)
adata_anno = st.io.stereo_to_anndata(data,flavor='seurat',output='out.h5ad')
