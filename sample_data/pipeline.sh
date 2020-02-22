
cd ..

# Create Segments and events-to-segments mapping
python yanagi.py segment -l 100 -wd sample_data/ -o segs_100 -ioe sample_data/AS_events.ioe


# Run yanagiCount
## Quasiindex
mkdir sample_data/quasiindex/
$RAPMAPPATH/rapmap quasiindex -i sample_data/quasiindex/ -t sample_data/segs_100.fa --segments sample_data/segs_100.fa.meta -k 23 --keepDuplicates

## Quasimap
mkdir sample_data/output/
$RAPMAPPATH/rapmap quasimap -i sample_data/quasiindex/ -s -1 sample_data/sample_1.fasta -2 sample_data/sample_2.fasta -o sample_data/output/ --consensusSlack 0.5 --noOrphans --hardFilter -t 10

# Calculate PSI values
mkdir sample_data/psis/
python yanagi.py psiCalc -es sample_data/segs_100_AS_events.ioe.flex.evs2segs -s sample_data/segs_100.fa.meta -i sample_data/output/ -o sample_data/psis/

# Run Tx Quantification
