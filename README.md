# CIN_HPV
Code accompanying the paper "HPV-induced host epigenetic reprogramming is lost upon progression to high-grade cervical intraepithelial neoplasia" accepted for publication in the *International Journal of Cancer*. <br>

## Scope and data availability
Two epigenome datasets were generated with Illumina MethlyationEPIC v1.0 bead chips. <br>
The first dataset was also used for the development of the WID-CIN score predicting cervical intraepithelial neoplasia grade 3 (see https://doi.org/10.1186/s13073-022-01116-9). <br>
The original raw microarray data is available in the European Genome-phenome Archive (EGA) database under accession code EGAS00001005078. <br>
The DNA used for the epigenome analysis was extracted from cervical smear samples from 5 subgroups of women: <br>
1. no cytological abnormalities, no high-risk HPV infection (HPV-)
2. no cytological abnormalities, high-risk HPV infection (HPV+)
3. cervical intraepithelial neoplasia grade 1 (CIN1)
4. cervical intraepithelial neoplasia grade 2 (CIN2)
5. cervical intraepithelial neoplasia grade 3 or invasive cervical cancer (CIN3+)

The second dataset was generated specifically for this study with DNA from an experiment with cancer cell lines where apoptosis was induced using BH3 mimetics. <br>
The original raw microarray data is available at the NCBI Gene Repository Omnibus (GEO) repository under accession code GSE207224. <br>

The analysis presented in the paper focusses on <br>
1) the genome-wide methylation changes in the host genome upon infection with high-risk HPV <br>
2) devising a classifier predicting hrHPV status, we term the WID-HPV <br>
3) investigating correlations of the WID-HPV with increasing CIN grade, replicative age and apoptosis.
