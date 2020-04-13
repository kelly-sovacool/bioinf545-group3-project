#!/bin/bash
wget ftp://ftp.cngb.org/pub/SciRAID/Microbiome/humanGut_9.9M/GeneCatalog/IGC.fa.gz
tar -xvf IGC.fa.gz
wget ftp://ftp.cngb.org/pub/SciRAID/Microbiome/humanGut_9.9M/GeneAnnotation/IGC.annotation_OF.summary.gz
gzip -d IGC.annotation_OF.summary.gz
