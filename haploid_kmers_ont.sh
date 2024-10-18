#!/bin/bash

H1=assembly_h1.fa
H2=assembly_h2.fa
K=17
ONT=ont.fq.gz

meryl k=$K count output h1.meryl $H1
meryl k=$K count output h2.meryl $H2
meryl difference h1.meryl h2.meryl output h1_only.meryl
meryl difference h2.meryl h1.meryl output h2_only.meryl
meryl-lookup -existence -sequence $ONT -mers h1_only.meryl h2_only.meryl -output ont_haploid_kmers.txt
