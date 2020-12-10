#!/usr/bin/env nextflow

/*
How to run:
nextflow -C phylogenetic_tree_nuc_seq.nf.config run phylogenetic_tree_nuc_seq.nf -profile amanj
*/

params.contigs='contig_files'
params.req='req_files'

/* input files */
//contig sequences
contig_files = Channel.fromFilePairs("${params.contigs}/*/*_filt_contigs.fa",size:1)

//tsv table 
tsv_file = Channel.fromFilePairs("${params.req}/*.tsv",size:1)

//nexus commands
nex_file = Channel.fromFilePairs("${params.req}/*.nex",size:1)

//sequences with outgroups
outgroup_file = Channel.fromFilePairs("${params.req}/*.fa",size:1)

/**
processes - contig data
**/

process extracting_results_based_on_keyword{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All", mode:'link'

  input:
  set name, tsv_table from tsv_file 
  
  output:
  set "tsv_1.txt" into tsv_filtering
  
  script:
""" 
 cat ${tsv_table[0]} | grep -i "anello\\|TTV\\|torque" | grep -iv "simian\\|mosq\\|chimp\\|gorilla\\|rodent\\|panda\\|paguma\\|seal\\|pine\\|porcine\\|tick" > "tsv_1.txt"
"""
}

combine_tsv_with_seq = contig_files.combine(tsv_filtering)

process selecting_seq_based_on_tsv{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/selecting_seq_based_on_tsv", mode:'link'

  input:
  set sample_id, sequences, tsv_table from combine_tsv_with_seq 
  
  output:
  set sample_id, "${sample_id}_list_with_selected_contigs.txt", "${sample_id}_extracted_contigs.fa" into extracted_contigs
  
  script:
""" 
 cat ${tsv_table} | grep '${sample_id}' | awk '{print \$2}' > "${sample_id}_list_with_selected_contigs.txt"
 seqtk subseq  ${sequences[0]}  ${sample_id}_list_with_selected_contigs.txt  > "${sample_id}_extracted_contigs.fa"
"""

}

process filter_contigs{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/selecting_seq_based_on_tsv", mode:'link'

  input:
  set sample_id, list, contig from extracted_contigs 
  
  output:
  set sample_id, stdout, "${sample_id}_selected_seq_filt.fa" into filt_contigs
  
  script:
""" 
 seqtk seq -L ${params.seq_size} ${contig} > "${sample_id}_selected_seq_filt.fa"
  nr_of_seq=\$(cat "${sample_id}_selected_seq_filt.fa" | grep ">" | wc -l)
 if [ \$nr_of_seq -eq 0 ]
 then
    echo false
 else
    echo true
 fi
"""
}

combine_filt_contigs_outgroup = filt_contigs.filter({ it[1].contains("true") }).combine(outgroup_file)

process add_outgroup_to_seq{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/outgroup_added", mode:'link'

  input:
  set sample_id, bool, contig, outgrp_name, outgrp_seq from combine_filt_contigs_outgroup 
  
  output:
  set sample_id,"${sample_id}_selected_seq_filt_with_outgrp.fa" into filt_contigs_and_outgrp
  
  script:
""" 
 cat ${contig} ${outgrp_seq[0]} > "${sample_id}_selected_seq_filt_with_outgrp.fa"
"""
}

process mafft{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/mafft", mode:'link'

  input:
  set sample_id, seq from filt_contigs_and_outgrp 
  
  output:
  set sample_id,"${sample_id}_mafft_alignment.fasta" into mafft_out
  
  script:
""" 
 mafft --thread ${task.cpus} --threadtb 5 --threadit 0 --inputorder --anysymbol --auto ${seq} > "${sample_id}_mafft_alignment.fasta" 
"""
}

process trimAl{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/trimAl", mode:'link'

  input:
  set sample_id, mafft_aln from mafft_out 
  
  output:
  set sample_id, "${sample_id}_trimAl_trimmed.fasta" into trim_out
  
  script:
""" 
 trimal -in ${mafft_aln} -out "${sample_id}_trimAl_trimmed.fasta" -automated1
"""
}

trim_out.into{raxml_in;nexus_in}
combine_nexus_commands = nexus_in.combine(nex_file)


process create_nexus_file{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/nexus", mode:'link'

  input:
  set sample_id, seq, command_file_name, nex_file from combine_nexus_commands 
  
  output:
  set sample_id, "${sample_id}_.nex" into nex_out
  
  script:
""" 
 seqmagick convert --output-format nexus --alphabet dna ${seq} ${sample_id}_.nex
 cat ${nex_file[0]} >> ${sample_id}_.nex
"""
}

process RAxML{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/RAxML", mode:'link'

  input:
  set sample_id, seq from raxml_in 
  
  output:
  set sample_id, "${sample_id}*" into raxml_out
  
  script:
""" 
  raxml-ng --msa ${seq} --model GTR+FO+G --opt-branches on --opt-model on --tree pars{10},rand{10} --all --bs-trees 100 --force --threads ${task.cpus} --prefix ${sample_id}

"""
}
/*
process mrbayes{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/mrbayes", mode:'link'

  input:
  set sample_id, seq_nex from nex_out 
  
  output:
  set sample_id, "${sample_id}*", 'log.txt' into mrbayes_out
  
  script:
""" 
  cat ${seq_nex} > "${sample_id}.nex"
  mpirun -np ${task.cpus} mb -i "${sample_id}.nex"
"""
}
*/
