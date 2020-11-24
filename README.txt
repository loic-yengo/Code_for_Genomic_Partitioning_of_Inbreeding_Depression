We provide a general overview of our analysis pipeline.
The scripts provided here are now readily useable.
However, readers/users could easily borrow pieces of this code to run their own analyses.

1) Calculate annotation-specific inbreeding measures

For FUNI, that step utilises multiple scripts as described below
a- qc_annot.R which generates qced annotations
b- fannot (fannot.cpp + Makefile) which is called from the pbs script: runFannot_per_block.pbs

For FROJH, we for annotate each of the ROH detected (prepare_input_for_XROH.R, XROH.R)
then merge the results using combine_ROH.R. Data for our sensitivity analysis using increasing length of ROH
are generated using the script combine_ROH_sensitivity.R.

2) Use summary statistics to estimate enrichment of heritability and that of inbreeding depression (ID).
Once the GWAS are performed we prepare the data for analyses using prepare_gwas_data.R, then run the actual analysis using enrichment_h2_ID.R.
Results, for this analyses correspond to Table S6.

3) Estimation of ID using LD score regression is described in ID-LDSC.R.

4) Example_script_for_Estimation_of_ID_enrichment.R gives a detailed example of how run our model to estimate ID enrichment. Try it!
