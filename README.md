# ST_upland_rice
Scripts used for spatial transcriptomics of upland rice. 

1.The modified gff file was generated using gtf.generator in gff_edit/

2.Data integration was performed using merge.sh, and RDS files generated using Harmony were used in downstream analysis.

3.The marker genes were found using CR_cal_marker/cal_marker.sh

4. Following analysis were performed using an adjusted RDS which was generated using diff_gene_SCT_position_adjusting/position_adjusting.R

5.The RNA velocity analysis was performed by using scvelo/merge_analysis.py

6.Pseudotime trajectory analysis were performed by using CR_monocle/monocle3.R, while positional trajectory analysis were performed by using CR_monocle/monocle3.2.R (coleoptilar node dataset) or monocle_position/monocle_position.R (root tip dataset). crown root primodrdia were selected by using CR_monocle/S00_CR_selected.R

7. Genes expressed in same vertical region were find by using vertical_cor_genes/work.sh

8.DEGs between upland rice and irrigated were find by using diff_gene_SCT_position_adjusting/diff_gene.sh

9.The exhibition of hormone related genes were generated by using known_marker.R

