# Analytic code: Abstraction of Longitudinal Breast Oncology Records Using a Retrieval Pipeline with General-Purpose Large Language Models

R code that produces the tables, figures and statistics reported in the manuscript (under review at npj Digital Medicine).
The patient-level data contain protected health information and are not included. Before running, set the
`JUNIOR_FIRST_PAPER_BASE_DIR` environment variable to the study data folder; `config.R` reads it.

Run order:
1. `cleanup/` scripts: clean the expert, research coordinator and LLM outputs into the intermediate files.
2. `NonDrugAnalysis.R` and `DrugComparison.R`: per-variable and per-drug comparisons.
3. `manuscript_tables.R` (sources `line_metrics.R`) and `manuscript_figures.R`: manuscript and supplementary tables and figures.
4. `paired_bootstrap_differences.R .`: paired patient-level bootstrap differences between sources: each LLM vs the
   second oncologist, the two best LLMs vs the two oncologists, research coordinators vs each LLM, median survival,
   and the ratios of hazard ratios in Supplementary Table 9.

The pipeline that produced the LLM abstractions is at https://github.com/MrJimb0/junior
