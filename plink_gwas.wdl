version 1.0

## Version 07-08-2021
##
## This workflow runs a PLINK GWAS.
## PLINK2 documentation: https://www.cog-genomics.org/plink/2.0/
##
## Cromwell version support - Successfully tested on v65
##
## Distributed under terms of the MIT License
## Copyright (c) 2021 Brian Sharber
## Contact <brian.sharber@vumc.org>

workflow plink_gwas {
    input {
        File input_parameters_file    # File containing tab-separated list of: bed file, bim file, bam file, and output name for the file at the end of the analysis.
        String plink2_docker = "briansha/plink2:terra"
        String r_base_docker = "briansha/plink_gwas:latest"
    }
    Array[Array[String]] input_parameters = read_tsv(input_parameters_file)

    scatter (line in input_parameters) {
        call GWAS {
          input:
            bed_file = line[0],
            bim_file = line[1],
            fam_file = line[2],
            output_name = line[3],
            docker = plink2_docker
        }
    }
    call GWASAnalysis {
      input:
        input_files = GWAS.output_file, # Refers implicitly to the entire array of output files from GWAS
        docker = r_base_docker
    }

    output {
        Array[File] output_file = GWAS.output_file
        File ADD_file = GWASAnalysis.ADD_file
        File ADD_MCA_file = GWASAnalysis.ADD_MCA_file
        Array[File] plots = GWASAnalysis.plots
    }

    meta {
    	author : "Brian Sharber"
        email : "brian.sharber@vumc.org"
        description : "This workflow runs a PLINK GWAS."
    }
}

task GWAS {
    input {
        File bed_file
        File bim_file
        File fam_file
        String output_name

        # PLINK2 parameters
        File pheno
        String pheno_name
        File covar
        String covar_name
        String parameters

        # Runtime
        String docker # Docker image containing PLINK2.
        Int? disk_size_override
        Float memory = 12.0
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float bed_size = size(bed_file, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * bed_size)])
    Float memory_for_plink = ceil(memory * 1000)

    command <<<
        set -euo pipefail

        plink2 \
        --bed ~{bed_file} \
        --bim ~{bim_file} \
        --fam ~{fam_file} \
        --logistic interaction  \
        --out ~{output_name} \
        --pheno ~{pheno} \
        --pheno-name ~{pheno_name} \
        --1 \
        --covar ~{covar} \
        --covar-name ~{covar_name} \
        --covar-variance-standardize \
        --maf 0.05 \
        --threads ~{cpu} \
        --memory ~{memory_for_plink} \
        --parameters ~{parameters}
    >>>

    output {
        File output_file = "${output_name}.${pheno_name}.glm.logistic"
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}

task GWASAnalysis {
    input {
        Array[File] input_files
        String add_colname = "ADD" # Column name for the ADD data
        String addx_colname        # Column name for the ADDx data

        # Runtime
        String docker # Docker image containing PLINK2.
        Int? disk_size_override
        Float memory = 12.0
        Int cpu = 1
        Int preemptible = 1
        Int maxRetries = 0
    }
    Float input_size = size(input_files, "GiB")
    Int disk = select_first([disk_size_override, ceil(10.0 + 2.0 * input_size)])

    command <<<
        set -euo pipefail

        R --no-save --args ~{add_colname} ~{addx_colname} ~{sep=' ' input_files} <<RSCRIPT
        library(data.table)
        args <- commandArgs(trailingOnly = TRUE)
        chr <- 1
        count <- 0
        # Skip the first two args - they are column indicators, not files.
        for (file in args) {
          if (count < 2) {
            count <- count + 1
            }
            else {
              gwas <-fread(file)
              subset <-subset.data.frame(gwas, TEST==args[1])
              complete_cases <-subset[complete.cases(subset), ]
              colnames(complete_cases)[1]<-"CHROM"
              complete_cases <-subset[complete.cases(subset), ]
              colnames(complete_cases)[1]<-"CHROM"
              subset1<-complete_cases[,c("CHROM","POS", "ID", "P")]
              subset1[,"P"] <- as.numeric(unlist(subset1[,"P"]))
              subset1[,"CHROM"] <-as.numeric(unlist(subset1[,"CHROM"]))
              subset1[,"POS"] <-as.numeric(unlist(subset1[,"POS"]))
              subset22 <-subset1[complete.cases(subset1), ]
              output_file <- paste(chr, "ADD.txt", sep="_")
              write.table(subset22, output_file, col.names=T, row.names=F, quote=F, sep="\t")

              subset <-subset.data.frame(gwas, TEST==args[2])
              complete_cases <-subset[complete.cases(subset), ]
              colnames(complete_cases)[1]<-"CHROM"
              complete_cases <-subset[complete.cases(subset), ]
              colnames(complete_cases)[1]<-"CHROM"
              subset1<-complete_cases[,c("CHROM","POS", "ID", "P")]
              subset1[,"P"] <- as.numeric(unlist(subset1[,"P"]))
              subset1[,"CHROM"] <-as.numeric(unlist(subset1[,"CHROM"]))
              subset1[,"POS"] <-as.numeric(unlist(subset1[,"POS"]))
              subset22 <-subset1[complete.cases(subset1), ]
              output_file <- paste(chr, "ADDxnew_MCA.txt", sep="_")
              write.table(subset22, output_file, col.names=T, row.names=F, quote=F, sep="\t")
              chr <- chr + 1
            }
        }
        RSCRIPT

        cat *ADD.txt >> ADD_all_chr_results.txt
        cat *ADDxnew_MCA.txt >> ~{addx_colname}_all_chr_results.txt

        R --no-save --args ~{addx_colname} ~{addx_colname}_all_chr_results.txt <<RSCRIPT
        library(data.table)
        library(qqman)
        args <- commandArgs(trailingOnly = TRUE)
        count <- 0
        for (file in args) {
          if (count < 1) {
            count <- count + 1
            }
          else {
            ADD <-fread("ADD_all_chr_results.txt")
            ADDxMCA <-fread(file)
            ADD_subset <-subset.data.frame(ADD)
            ADD_subset[,"P"] <- as.numeric(unlist(ADD_subset[,"P"]))
            ADD_subset[,"CHROM"] <-as.numeric(unlist(ADD_subset[,"CHROM"]))
            ADD_subset[,"POS"] <-as.numeric(unlist(ADD_subset[,"POS"]))
            ADD_subset <- ADD_subset[complete.cases(ADD_subset), ]
            ADDxMCA_subset <-subset.data.frame(ADDxMCA)
            ADDxMCA_subset[,"P"] <- as.numeric(unlist(ADDxMCA_subset[,"P"]))
            ADDxMCA_subset[,"CHROM"] <-as.numeric(unlist(ADDxMCA_subset[,"CHROM"]))
            ADDxMCA_subset[,"POS"] <-as.numeric(unlist(ADDxMCA_subset[,"POS"]))
            ADDxMCA_subset <- ADDxMCA_subset[complete.cases(ADDxMCA_subset), ]
            png("ADD_all_chr_results_qqplot.png")
            print(qq(as.numeric(unlist(ADD[,"P"]))))
            dev.off()
            png("ADD_all_chr_results_manhattan.png")
            print(manhattan(ADD_subset, chr="CHROM", bp="POS", snp="ID", p="P", annotatePval = 1E-5, ylim=c(0,20)))
            dev.off()
            ADDx_qq_output_name <- paste(args[1], "_qq.png", sep="")
            ADDx_manhattan_output_name <- paste(args[1], "_manhattan.png", sep="")
            png(ADDx_qq_output_name)
            print(qq(as.numeric(unlist(ADDxMCA[,"P"]))))
            dev.off()
            png(ADDx_manhattan_output_name)
            print(manhattan(ADDxMCA_subset, chr="CHROM", bp="POS", snp="ID", p="P", annotatePval = 1E-5, ylim=c(0,20)))
            dev.off()
          }
        }
        RSCRIPT
    >>>

    output {
        File ADD_file = "ADD_all_chr_results.txt"
        File ADD_MCA_file = "${addx_colname}_all_chr_results.txt"
        Array[File] plots = glob("*.png")
    }

    runtime {
        docker: docker
        memory: memory + " GiB"
		disks: "local-disk " + disk + " HDD"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}
