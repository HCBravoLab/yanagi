#!/bin/bash -x

# Merge all .ioe files and remove duplicate headers
cat test_hg19_71/events/events_5types_*.ioe | awk '/seqname/&&c++>0 {next} 1' - > test_hg19_71/hg19_71_events.ioe

# Remove patch and decoy sequences, to only keep events on the primary assembly
sed -i '/PATCH/d;/HSCHR/d;/NOVEL/d;/GL/d' test_hg19_71/hg19_71_events.ioe