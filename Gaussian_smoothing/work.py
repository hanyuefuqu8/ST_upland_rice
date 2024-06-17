#the first parameter is a directory of gem files, the second parameter is genelist
import stereo as st
import os
from matplotlib.pyplot import savefig
import PIL.Image as Image
import matplotlib.pyplot as plt
import pandas as pd
import sys

def Gaussian(path,ID,name):
  toImage = Image.new('RGB', (2160 * 3, 640 * 10))
  count = 0
  begin_x = 0
  begin_y = 0
  dir=path
  files=os.listdir(dir)
  for file_name in files:
      if file_name.endswith('.gem'):
        fname = os.path.join(dir, file_name)
        data=st.io.read_gem(fname,bin_size=40)
        data.tl.cal_qc()
        data.tl.filter_cells(min_gene=0, pct_counts_mt=0)
        data.tl.filter_genes(min_cell=0)
        data.tl.raw_checkpoint()
        data.tl.normalize_total(target_sum=10000)
        data.tl.log1p()
        data.tl.pca(use_highly_genes=False, n_pcs=50, svd_solver='arpack')
        data.tl.gaussian_smooth(n_neighbors=10, smooth_threshold=90)
        data.tl.scale(max_value=10) #only for gaussian_smooth_scatter_by_gene
        try:
          data.plt.gaussian_smooth_scatter_by_gene(gene_name=ID,dot_size=30)
        except TypeError as e:
          print(e)
        except Exception as e:
          print(e)
        filename=name + file_name + '.pdf'
        savefig(filename)
       # img=Image.open(filename)
       # toImage.paste(img, (begin_x, begin_y))
       # begin_x += 2160
       # if begin_x % 6480 < 2160 and begin_x != 0:
        #      begin_x = 0
        #      begin_y += 640
#  filename=name + file_name + '.pdf'
#  toImage.save(filename)


if __name__ == '__main__':
    """main
    """
    data=pd.read_csv(sys.argv[2])
    path=sys.argv[1]
    for i in range(len(data)):
        Gaussian(path,data.iat[i,0],data.iat[i,1])
