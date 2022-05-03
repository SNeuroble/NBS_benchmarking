#!/bin/bash
#
# convert CompROI output into more simple readable format for loading into NBS script
# also renam accordingly
#
# existing format: project/HCP_S1200_activation/by_nets/SOCIAL_cope6/GroupSize484/lower_level/100307_cope6_by_net.csv
# target format:   project/HCP_S1200_activation/by_nets/SOCIAL/100307_SOCIAL_cope.txt
# template format:  project/HCP_S1200/SOCIAL/100206_SOCIAL_GSR_matrix.txt
#
# next steps:
#    change in setparams_bench: data_type_suffix
#    change in setpaths: input data-> by_nets

base_base_data_dir="/home/smn33/project/HCP_S1200_activation/by_nets/"
task="SOCIAL"
cope="6"

base_data_dir="$base_base_data_dir/$task"
in_data_dir="${base_data_dir}_cope${cope}/GroupSize484/lower_level"
out_data_dir="$base_data_dir"

subIDs_suffix='_subIDs.txt';
infile_suffix="_cope${cope}_by_net.csv"
outfile_suffix="_cope.txt"

subIDs_file="${out_data_dir}$subIDs_suffix"

tmp=$(mktemp)

mkdir -p $out_data_dir

# make subIDs file
ls $in_data_dir > $subIDs_file
sed '/^t/d' $subIDs_file > $tmp
cut -d "_" -f1 $tmp > $subIDs_file

# convert all sub data:
# remove first column, remove first row, commas->tabs
while read subject; do 
    infile="$in_data_dir/${subject}$infile_suffix"
    outfile="$out_data_dir/${subject}_${task}$outfile_suffix"
    cut -d "," -f3- $infile | sed -n '2p' | sed 's/,/\n/g' > $outfile
    #printf "$infile\n$outfile\n$subIDs_file\n"
done < $subIDs_file
    


