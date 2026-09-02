# Metatranscriptome processing pipeline - read/adaptor trimming, de novo assembly, annotation, and quantification

### Step 1: Clone this repository

```
git clone https://github.com/annikagomez/dyhrman_lab_metaT_pipeline.git
cd dyhrman_lab_metaT_pipeline
```

### Step 2: Install Snakemake using conda
- Create and activate conda environment:
```
conda create --name snakemake
conda activate snakemake
```

- Install Snakemake (I built this pipeline with v8.25.5):
```
conda install -c bioconda snakemake=8.25.5
```

### Step 3: Edit Snakefile, slurm script, and config.v8+.yaml
- Snakefile global wildcards: In line 1 of the Snakefile, provide the absolute path of the folder containing your raw sequences and make sure the file ending (_R1_001.fastq.gz) matches that of your forward read files (SeqCenter files usually have the _R1_001.fq.gz pattern). We only include the forward reads because we want a list of all of the sample names (1 for each sample), not all of the files (2 for each sample, forward and reverse reads)
- Snakefile trim_galore input: In the input to rule trim_galore, edit r1 and r2 to match your input location. r1 should match the global wildcard in line 1 and r2 should point to the reverse reads
  - To skip any of the steps, add the path to the reads you want to start with as global wildcards and edit the input of the rule you want to start with to be the global wildcards
- profile/config.v8+.yaml: This file tells snakemake what default slurm resources to use to run a job. These resources will be used unless you specificy different resource parameters in a rule. Add the account you want to pull the resources from. You can adjust any of the other parameters as well - eg I have it set to only submit 5 batch jobs at once.
- snakemake.sh: Add the account you want to pull resources from to run Snakemake (doesn't need much) and activate your snakemake environment

### Step 4: Unzip rRNA reference file
```
gunzip ref_seqs/smr_v4.3_default_db.fasta.gz
```

### Step 5: Build DAG graph
- This step visualizes the workflow Snakemake will run and is a good way to troubleshoot any issues without submitting the slurm script and waiting for it to run
```
snakemake --forceall --dag | dot -Tpdf > dag.pdf
```
- Look at the output PDF to make sure Snakemake will do what you expect

### Step 6: Submit slurm script to  run Snakemake
```
sbatch snakemake.sh
```
