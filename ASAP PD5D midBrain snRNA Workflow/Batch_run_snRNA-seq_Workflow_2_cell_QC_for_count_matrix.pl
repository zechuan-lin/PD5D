#!usr/bin/perl
#### Put this script and the R script Single_cell_QC_for_count_matrix.R under the folder storing CellRanger output folder of samples, change the file matching pattern below based on your file names ######################
@ARGV = <batch*\_BN*>;
my $i = 1;
foreach(@ARGV){
        system "Rscript Single_cell_QC_for_count_matrix.R $_";
}
