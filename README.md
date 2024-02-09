# AmplificationTimeR_paper_archive
A repository to store any code that was run as part of the `AmplificationTimeR` publication.

## Overview of files and folders
- `Generating_figures`: A folder containing Rmarkdown documents for generating figures that form part of the `AmplificationTimeR` manuscript, as well as their html output, and session info files listing tools and programme versions.

    | File | Description |
    | ---- | ----------- |
    |`AmplificationTimeR_paper_figures_tables_analysis_080224.Rmd` | Rmarkdown file with analysis and figures generated using PCAWG timing data. |
    | `AmplificationTimeR_paper_figures_tables_analysis_080224.html` | html file generated from the above Rmarkdown document. |
    | `AmplificationTimeR_plotting_PCAWG_data_sessionInfo_2024-02-09.txt` | session info for the above Rmarkdown and associated html file. |
    | `plotting_cov_c8_wgd_csf3.Rmd` | Rmarkdown document for analysis and plotting of 8q segments timed. |
    | `plotting_cov_c8_wgd_csf3.html` | html file generated from the above Rmarkdown document. |
    | `plotting_cov_c8_sessionInfo_2024-02-08.txt` | session info for the above Rmarkdown and associated html file. |
    | `plotting_output_of_simulations_05_08.01.24_csf3.Rmd` | Rmarkdown document for plotting figures from simulated data. |
    | `plotting_output_of_simulations_05_08.01.24_csf3.nb.html` | html file associated with the above Rmarkdown document. |
    | `multinomial_plotting_simulation_sessionInfo_2024-01-08.txt` | session info for the above Rmarkdown and html file. |

- `running_AmplificationTimeR`: A folder containing code associated with running `AmplificationTimeR` on a portion of the PCAWG data set.

    | File | Description |
    | ---- | ----------- |
    | `20231123_running_AmplificationTimeR_on_PCAWG.R` | R script used to run `AmplificationTimeR` on PCAWG data. |
    | `AmplificationTimeR_run_sessionInfo_2024-01-06.txt` | session info from most recent run of `AmplificationTimeR` on PCAWG data. |

- `running_cancerTiming`: A folder containing code associated with running `cancerTiming` on a portion of the PCAWG data set.

    | File | Description |
    | ---- | ----------- |
    | `20230109_cancerTiming_PCAWG_MYC_BRCA_OV_PRAD.R`| R script used to run `cancerTiming` on PCAWG data. |
    | `sessionInfo_2023-01-09.txt`| session info from the most recent run of `cancerTiming` on PCAWG data. |

- `running_MutationTimeR`: A folder containing code associated with running `MutationTimeR` on a portion of the PCAWG data set.

    | File | Description |
    | ---- | ----------- |
    | `20230110_MutationTimeR_PCAWG_BRCA_OV_PRAD_all_array.R` | R script used to run `MutationTimeR` on PCAWG data. |
    | `collating_MutationTimeR_results_250123.R` | R script to collate `MutationTimeR` timing data for segments spanning MYC. |
    | `sessionInfo_2023-01-19.txt` | session info from the most recent run of `MutationTimeR` on PCAWG data. |

- `Sample_mismatches`: 

    | File | Description |
    | ---- | ----------- |
    | `checking_sample_mismatches_070224.R` | R script used to for identify and check inconsistencies in data between Battenberg and DPClust files in PCAWG data. |

- `SBS_classification_of_mutations`:

    | File | Description |
    | ---- | ----------- |
    | `20221118_get_individual_mut_sigs_PCAWG.R` | R script for working out mutational signature of individual mutations. |

- `Testing_with_simulated_data`: A folder containing code associated with simulating data to test `AmplificationTimeR`

    | File | Description |
    | ---- | ----------- |
    | `testing_AmplificationTimeR_synthetic_multinomial_sampled_data_csf3.Rmd` | Rmarkdown document for simulating timing data and testing the effect of varying mutation number and clocklike proportion. |
    | `testing_AmplificationTimeR_synthetic_multinomial_sampled_data_csf3.html` | html file generated from the above markdown document. |
    | `multinomial_sampled_data_sessionInfo2024-01-08.txt` | session info for the Rmarkdown document and html file above. |
    | `testing_AmplificationTimeR_synthetic_mismatched_equations.Rmd` | Rmarkdown document for simulating timing data and testing effect of applying mismatched equations. |
    | `testing_AmplificationTimeR_synthetic_mismatched_equations.nb.html` | html file generated from the above Rmarkdown document. |
    | `multinomial_data_right_and_wrong_equations_sessionInfo_2024-01-05.txt` | session info file for the above Rmarkdown and html document. |

- `Timing_8q_segments`: A folder containing code for identifying and timing 8q segments for a portion of the PCAWG data set.

    | File | Description |
    | ---- | ----------- |
    | `20231116_list_all_8q_segments_BRCA_OV_PRAD.R` | An Rscript to list all segments on 8q in the PCAWG BRCA, OV, and PRAD data. |
    | `20231123_time_all_segments_on_8q_sbs1_sbs5_BRCA_OV_PRAD.R` | R script to time all 8q segments in the aforementioned PCAWG data. |
    | `BRCA_OV_PRAD_AmplificationTimeR_C8_results_sbs1_sbs5_sessionInfo2024-01-05.txt` | session info associated with timing of the aforementioned 8q segments. |
