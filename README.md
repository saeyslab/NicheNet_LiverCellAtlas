<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- github markdown built using
rmarkdown::render("README.Rmd",output_format = "md_document")
-->

# Differential NicheNet code from Liver Atlas paper (Guilliams et al., Cell 2022)

In this repository you can find all the code that was used for the
Differential NicheNet analyses and visualizations as described in the
Liver Atlas paper from Guilliams et al.: [Spatial proteogenomics reveals
distinct and evolutionarily conserved hepatic macrophage
niches](https://www.sciencedirect.com/science/article/pii/S0092867421014811).

**Differential NicheNet** is an extension of the default
[NicheNet](https://github.com/saeyslab/nichenetr) pipeline to compare
cell-cell interactions between different niches and better predict
niche-specific ligand-receptor (L-R) pairs. It was used in that paper to
predict ligand-receptor pairs specific for the Kupffer cell niche in
mouse and human.

For users interested in applying Differential NicheNet to their own
data, we recommend the vignette on the Github page of the nichenetr
package [Differential NicheNet analysis between niches of
interest](https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md)
or [Differential NicheNet analysis between conditions of
interest](https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet_pEMT.md).

## Content of this repository

-   scripts 0 and 1: scripts used to make Seurat objects from raw counts
    and cell annotations for the mouse data.
-   scripts 2a, b and c: scripts used to annotate mouse stellate cells,
    LSECs and hepatocytes as portal vein or central vein cells.
-   scripts 3a, b: scripts used to add this spatial zonation information
    to the existing Seurat object.
-   scripts 4a, b, and c: scripts used to run the Differential NicheNet
    pipeline on the mouse data, comparing the Kupffer cell niche
    separately to respectively the bile duct macrophage niche, capsule
    macrophage niche and central vein macrophage niche.
-   script 4d: script to analyze and combine these previous Differential
    NicheNet analyses, to come to a list of ligand-receptor pairs
    specific to the Kupffer cell niche (including code to make Figure
    S8G).
-   script 5: script used to run the Differential NicheNet pipeline on
    the mouse data, comparing Kupffer cell niche to the niches of all
    the other liver macrophage niches combined.
-   scripts 6 and 7: scripts used to make Seurat objects from raw counts
    and cell annotations for the human data.
-   script 8: script used to run the Differential NicheNet pipeline on
    the human data, comparing Kupffer cell niche to the niches of all
    the other liver macrophage niches combined.
-   script 9: script used to compare the mouse and human Differential
    NicheNet results, and look at evolutionary conservation, based on
    the output of scripts 5 and 8 (including code to make Figure 7B).
-   script 10: script used to predict upstream ligands driving
    differences in KC-specific genes between human and mouse.
-   script 11: script used to define and visualize the human-KC-specific
    ligand-receptor pairs, based on the output of scripts 9 and 10
    (including code to make Figure 6E).

Notes:

-   the exact data files you need to have in the data folder for running
    these scripts are available on this Zenodo page:
    <https://zenodo.org/record/5840333>. For exploration and downloading
    of other data from the paper, we refer to: [Liver Atlas Data
    Portal](https://www.livercellatlas.org/)

-   the code provided in this repository is tailored to work on our
    gridengine cluster (via the qsub package) and will need to be
    adapted to work on your system/server/cluster.

## References

Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular
communication by linking ligands to target genes. Nat Methods (2019)
<doi:10.1038/s41592-019-0667-5>

Guilliams et al.??Spatial proteogenomics reveals distinct and
evolutionarily conserved hepatic macrophage niches. Cell (2022)
<doi:10.1016/j.cell.2021.12.018>
