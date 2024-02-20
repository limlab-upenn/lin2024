# lin2024
codes used to analyze the impact of Kr dosage modulation on the transcriptional dynamics of eve

This depository contains Matlab files used in the manuscript "Multifaceted effects on even-skipped transcriptional dynamics upon Krüppel dosage changes"

eve_analysis
This file contains analysis of eve dynamics in each stripes and generate all the plots except false-color images (see Falsecolor) and plots regarding stripe position and stripe width (see eve_track_stripe). 

eve_track_stripe
This file track the stripe and width of individual eve stripes and genereate corresponding plots.

Falsecolor
This file generate images of embryos with active nuclei false-colored.

Fish_analysis
This file runs all the analysis and generate all plots for Fluorescence in situ hybridization (FISH) images.

Gap_analysis
This file contains analysis of gt dynamics and generate corresponding plots. Analysis of kni dynamics utilize similar code.

Gap_scaled_snapshot
This file generate snapshots of gt embryos with egg length scale and bars indicating expression domains.

Nuclei_width_embryo_length
This file measures the anterior-posterior (AP) and dorsal-ventral (DV) lengths of embryos, as well as measures the width of nucleus in turns of egg length (AP length).

Output_heatmap
Thia file generate heatmaps showing cumulative mRNA production in individual nuclei.

Shaded error bars were generated using code modified from: Rob Campbell (2023). raacampbell/shadedErrorBar (https://github.com/raacampbell/shadedErrorBar)


