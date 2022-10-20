threads=10
faa=../data/proteins/extracted_spikes/sarbecovirus_spikes.faa
#faa=../data/proteins/extracted_spikes/mers_spikes.faa

alignment_out=$(echo $faa|sed "s|.faa|.faa.aln|g"|sed "s|proteins/extracted_spikes|alignments/spike_protein_alignments|g")

#echo $faa
echo $alignment_out

mafft-linsi $faa > $alignment_out
