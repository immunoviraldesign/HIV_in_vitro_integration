# HIV_in_vitro_integration

During the first half of the viral life cycle, HIV-1 reverse transcribes its RNA genome and integrates the double-stranded DNA copy into a host cell chromosome. Despite progress in characterizing and inhibiting these processes, in situ mechanistic and structural studies remain challenging. This is because these operations are executed by individual viral preintegration complexes deep within cells. We therefore reconstituted and imaged the early stages of HIV-1 replication in a cell-free system. HIV-1 cores released from permeabilized virions supported efficient, capsid-dependent endogenous reverse transcription to produce double-stranded DNA genomes, which looped out from ruptured capsid walls. Concerted integration of both viral DNA ends into a target plasmid then proceeded in a cell extract-dependent reaction. This reconstituted system uncovers the role of the capsid in templating replication.

This GitHub repository contains the Python code used for identifying forward- and reverse-oriented integrations of HIV-1 into a target plasmid (pK184) and the R code used for plotting integration counts.  

Also included in this repository are the reference fasta sequences used for alignment:
verified_pk184.fa
brian_shifted_verifed_corrected_pk184.fa
HIV1_reference_crop_1000_9000.fa
HIV1_reference_5prime_LTR_to_N.fa

Processed data representing integration counts at each base in pK184 for three biological replicates:
046_insertions_directional_scaled_complex_CIGAR_final.xlsx
047_insertions_directional_scaled_complex_CIGAR_final.xlsx
048_insertions_directional_scaled_complex_CIGAR_final.xlsx

Processed data for ln-transformed integration counts for all three sequencing runs:
all_samples_insertions_directional_scaled_complex_CIGAR.xlsx

Integration sites mapped by Sanger sequencing of individual concerted integration clones:
concensus_integration_sites.txt

Bam files that contain 3’ integration site deep sequencing reads have been deposited in the NCBI Sequence Read Archive (SRA) database under Accession: PRJNA649355.

A description of the methods used to generate these data can be found in the forthcoming publication...

"Reconstitution and visualization of HIV-1 capsid-dependent replication and integration in vitro"

Authors:
Devin E. Christensen, Barbie K. Ganser-Pornillos, Jarrod S. Johnson, Owen Pornillos, Wesley I. Sundquist
