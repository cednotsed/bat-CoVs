threads=10
faa=../data/proteins/extracted_proteins/sarbecovirus_spikes.faa
faa=../data/proteins/extracted_proteins/MERS_spikes.faa

alignment_out=$(echo $faa|sed "s|.faa|.faa.aln|g"|sed "s|proteins/extracted_proteins|alignments|g")

#echo $faa
echo $alignment_out

mafft-linsi $faa > $alignment_out
