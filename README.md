# Proteomics-Differential-Expression-with-Fisher-exact-test
R script for Proteomics Differential Expression Analysis using Fisher exact test

## Analysis
* define (all possible) pairwise comparisons                                        
* filter proteins with low/high spectral counts                                     
* Fisher exact test for each protein                                                                
* p-value and multiple testing correction for each protein
* summarize results


## Input
* table (TAB delimited) as .txt, .csv, .xlsx (workbook with multiple sheets) 
* rows=spectral counts                                                             
* columns=conditions
  
## Output
* Excel workbook with one sheet for each pairwise comparison                          
* sheet content: Protein Identifier, Condition A, Condition B, log2 Fold Change, p-value, multiple testing adjusted
* p-value (annotation, etc.)
