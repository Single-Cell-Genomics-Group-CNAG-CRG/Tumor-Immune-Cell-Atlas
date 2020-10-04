#! /bin/bash

pipe_dir='/scratch/devel/melosua/phd/projects/TICA/misc/jobs/pipeline_dir'

pipe_fn="$pipe_dir/pipe_spotlight_oroph.pipe"
rm "$pipe_fn"
touch "$pipe_fn"

trn='melanoma'
for spatial in 161429 161430 161431 161432
do
printf "$spatial.decon\t$spatial\t-\t.\tn\t8:00:00\t1\t12\t.\tmodule purge; module load gsl/1.9_64 gsl/2.4 gcc/6.3.0 gmp/6.1.2 R/3.6.0 hdf5/1.10.1; analysis/aussie_oro/spotlight_deconv_oroph_job.R $trn $spatial\n" >> "$pipe_fn"
done
/home/devel/melosua/bin/cnag_pipeline.pl "$pipe_fn"
