#!/bin/bash

echo "[clean_up.sh Cleaning up temporary files...]
";

for file in RData phased_seg_sites.pos stitch.Herato1603.3450000.3550000.header stitch.Herato1603.3450000.3550000.PL.vcf.gz stitch.Herato1603.3450000.3550000.PL.vcf.gz.tbi focal.po focal.pos.match focal_TRUTH.pos_idx_correction tmp; do 
        rm -fr $file;
done
