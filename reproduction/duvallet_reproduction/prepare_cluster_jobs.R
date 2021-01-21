# ##############################################################################
#
# Function to prepare the cluster jobs for the reproduction of Duvallet et al.
#
# ##############################################################################

save.slurm <- function(data.tag, target=NULL, ...){
  
  # get other args, like time or memory (potential for optimization)
  other.args <- list(...)
  
  time <- ifelse('time' %in% names(other.args), other.args$time, '30:00')
  mem <- ifelse('mem' %in% names(other.args), other.args$mem, '5G')
  
  # if target is null, get target out of the data.tag
  if (is.null(target)){
    target <- strsplit(data.tag, split='_')[[1]][1]
  }
  
  # header for all
  slurm.script <- c("#!/bin/bash", 
                    paste0("#SBATCH -A zeller\n#SBATCH -t ", time), 
                    paste0("#SBATCH --mem ",  mem), 
                    paste0("#SBATCH -o /home/jawirbel/duvallet_reprod/results/", data.tag, '_', target, ".out"),
                    paste0("#SBATCH -e /home/jawirbel/duvallet_reprod/results/", data.tag, '_', target, ".err"),
                    "",
                    "module load  R/3.5.0-foss-2017b-X11-20171023",
                    "",
                    paste0("cp /home/jawirbel/microbiomeHD/data/clean_tables/", data.tag, '.metadata.clean.feather $TMPDIR/'),
                    paste0("cp /home/jawirbel/microbiomeHD/data/clean_tables/", data.tag, '.otu_table.clean.feather $TMPDIR/'),
                    "",
                    paste0("Rscript /home/jawirbel/duvallet_reprod/src/duvellet_reprod.R --metadata_in $TMPDIR/", data.tag, 
                           '.metadata.clean.feather --feat_in $TMPDIR/', data.tag, '.otu_table.clean.feather --label_in ', target
                    ),
                    "")
  file.out <- paste0('./execution/', data.tag, '_', target, '_slurm.sh')
  writeLines(slurm.script, con=file.out)
  
}



