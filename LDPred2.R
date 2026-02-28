library(bigsnpr)

rds_path <- snp_readBed("clean_plink_data.bed", backingfile="clean_plink_data")
obj.bigSNP <- snp_attach(rds_path)