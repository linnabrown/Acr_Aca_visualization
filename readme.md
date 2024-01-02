# Visualization Implementation for Acr-Aca interactions on gene circle and Crispr-Cas chain

Implemented in AcrDB website (https://bcb.unl.edu/AcrDB/) 

Example Page (https://bcb.unl.edu/AcrDB/anti_crispr_results.php?type=ncbi&organism=GCF_000015125.1)



## Visuzaltion Example

Circular genome view

1. The start/end position of the genome. The start and end point of the chromosome/contig is always in the middle lowest position of the circle. The end position is the total length of the chromosome/contig.

2. Acr-Aca operon/loci position. (Arrow Shape, "+" in brown, "-" in blue). Acr-Aca loci are visualized as arrows in brown/blue on the circle.

3. Spacer position. ( Ellipses in stroke purple color filled with white color). Only self-targeting spacers (STSs) are shown, i.e., a genome might contain multiple CRISPR-Cas loci, and those do not contain STSs are not shown here. Spacer position in CRISPR array is predicted by CRISPRCasFinder. STSs are shown as ellipses in stroke purple color filled with white color. The size is not proportionaly to the real spacer length. The labels indicating the actual postion of the STS (e.g. 1 Spacer Postion:976812-976847) are connected with a straight line to the ellipses on the circle.

4. Cas loci. ( Pink filled square with pink lables ). The Cas loci are predicted by CRISPRCasFinder. A genome might contain multiple CRISPR-Cas loci, and those do not contain STSs are not shown here.

5. Target postion. ( Ellipse in stroke in yellow green color filled with white color). The STS targets are identified by a BLASTn search with all the CRISPR spacers as query (see abovee in 1.1 Tabular view for GBA: column 20-21). The labels indicating the actual position of the STS targets (e.g. 2 Target Position:1946902-1946867) are connected with a straight line to the ellipses on the circle.

6. STS and target connection arcs. (The blue arc connecting spacers and targets inside the circle). We represent the STS spacer-target pairs by using arc lines in blue to connect them.
7. 
![image](https://github.com/linnabrown/Acr_Aca_visualization/assets/9990190/05152594-2fa6-46c5-9617-a202cdd22764)


## Citation

If you want to use this code for academic research, please cite us:

Le Huang, Bowen Yang, Haidong Yi, et al. AcrDB: a database of anti-CRISPR operons in prokaryotes and viruses. Nucleic Acids Res. 2021;49(D1):D622-D629. doi:10.1093/nar/gkaa857

```
@article{huang2021acrdb,
  title={AcrDB: a database of anti-CRISPR operons in prokaryotes and viruses},
  author={Huang, Le and Yang, Bowen and Yi, Haidong and Asif, Amina and Wang, Jiawei and Lithgow, Trevor and Zhang, Han and Minhas, Fayyaz ul Amir Afsar and Yin, Yanbin},
  journal={Nucleic Acids Research},
  volume={49},
  number={D1},
  pages={D622--D629},
  year={2021},
  publisher={Oxford University Press}
}
```
