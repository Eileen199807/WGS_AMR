#!/bin/bash
#####################################################################################
###                                                                                 #
### slurm-mpi_omp.cmd :                                                             #
### A SLURM submission script for running Hybrid MPI+OpenMP jobs on HPC2021 system  #
###                                                                                 #
### Compilation for Hybrid MPI+OpenMP program                                       #
### (1) Intel MPI Libaries                                                          #
###     module load impi                                                            #
###     mpiifort -qopenmp mpi_omp_hello.f90 -o mpi_omp_hello-f90-impi               #
### (2) OpenMPI Libaries                                                            #
###     module load openmpi                                                         #
###     mpif90 -fopenmp mpi_omp_hello.f90 -o mpi_omp_hello-f90-ompi                 #
###                                                                                 #
### Job submission                                                                  #
###    cd to directory containing program/executable, then:                         #
###    sbatch <location of this script>/slurm-mpi_omp.cmd                           #
###                                                                                 #
### SLURM for Multicore/Multi-thread: https://slurm.schedmd.com/mc_support.html     #
###                                                                                 #
### - Written by Lilian Chan, HKU ITS (2021-3-2)                                    #
###                                                                                 #
#####################################################################################
#SBATCH --job-name=trycycler
#SBATCH --mail-type=FAIL,END                      # Mail events
#SBATCH --mail-user=lwleen@hku.hk                 # Update your email address
#SBATCH --partition=intel                         # Specific Partition (intel/amd)
#SBATCH --time=5-10:00:00                         # Wall time limit (days-hrs:min:sec)
#SBATCH --ntasks=6                                # Total number of MPI tasks(processes)
#SBATCH --nodes=6                                 # Total number of compute node(s)
#SBATCH --ntasks-per-node=1                       # Number of MPI Tasks on each node
#SBATCH --cpus-per-task=4                         # Number of CPUs per each MPI task
#SBATCH --mem-per-cpu=6G                          # Memory setting
#SBATCH --output=%x.out.%j                        # Standard output file
#SBATCH --error=%x.err.%j                         # Standard error file
#SBATCH --account=sph_pengwu                      # Which group to use sbs_ssin or sph_pengwu
#####################################################################################
### The following stuff will be executed in the first allocated node.               #
### Please don't modify it                                                          #
#####################################################################################
echo "SLURMD_NODENAME       : $SLURMD_NODENAME"
echo "SLURM_NTASKS          : $SLURM_NTASKS"
echo "SLURM_JOB_NUM_NODES   : $SLURM_JOB_NUM_NODES"
echo "SLURM_CPUS_PER_TASK   : $SLURM_CPUS_PER_TASK"
echo "SLURM_CPUS_ON_NODE    : $SLURM_CPUS_ON_NODE"

echo JOBID ${SLURM_JOBID} : ${SLURM_NTASKS} CPUs allocated from ${SLURM_JOB_NODELIST}
echo Working directory is ${SLURM_SUBMIT_DIR}  1>&2
echo This SLURM script is running on host ${SLURMD_NODENAME} 1>&2

#################################################################################

##conda activate wgs-genetics


threads=32
nextdenovo_dir="/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/08_BENAGEN/02_ONT/02_assembly/02_devovo"
NextDenovo_dir="/lustre1/g/sph_pengwu/applications/NextDenovo"
NextPolish_dir="/lustre1/g/sph_pengwu/applications/NextPolish"
fq_dir="/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/08_BENAGEN/02_ONT/01_clean_reads"
subsample_dir="/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/08_BENAGEN/02_ONT/02_assembly/01_subsample"
minipolish_dir="/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/08_BENAGEN/CODE"
netcat_dir="/lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/Softwares/NECAT/Linux-amd64/bin"
genome_size="5500000"

for ID in `cat ECOLI.lst`
do  {
    mkdir ${subsample_dir}/${ID} 
    mkdir ${nextdenovo_dir}/${ID} 
    trycycler subsample --reads ${fq_dir}/${ID}.fastq.gz --out_dir ${subsample_dir}/${ID} --count 24 --genome_size "$genome_size"

    for i in 01 07 13 19; do
        canu -p canu -d canu_temp -fast genomeSize="$genome_size" useGrid=false -nanopore ${subsample_dir}/${ID}/sample_"$i".fastq
        python /lustre1/g/sph_pengwu/tutorial/2023_RIF_WGS_workshop/workshop3/Eileen/08_BENAGEN/CODE/canu_trim.py canu_temp/canu.contigs.fasta > ${nextdenovo_dir}/${ID}/assembly_"$i".fasta
        rm -rf canu_temp
    done

    for i in 02 08 14 20; do
        flye --nano-hq ${subsample_dir}/${ID}/sample_"$i".fastq --threads "$threads" --out-dir flye_temp
        cp flye_temp/assembly.fasta ${nextdenovo_dir}/${ID}/assembly_"$i".fasta
        cp flye_temp/assembly_graph.gfa ${nextdenovo_dir}/${ID}/assembly_"$i".gfa
        rm -rf flye_temp
    done

    for i in 03 09 15 21; do
        ${minipolish_dir}/miniasm_and_minipolish.sh ${subsample_dir}/${ID}/sample_"$i".fastq "$threads" > ${nextdenovo_dir}/${ID}/assembly_"$i".gfa
        any2fasta ${nextdenovo_dir}/${ID}/assembly_"$i".gfa > ${nextdenovo_dir}/${ID}/assembly_"$i".fasta
    done

    for i in 04 10 16 22; do
        ${netcat_dir}/necat.pl config config.txt
        realpath ${subsample_dir}/${ID}/sample_"$i".fastq > read_list.txt
        sed -i "s/PROJECT=/PROJECT=necat/" config.txt
        sed -i "s/ONT_READ_LIST=/ONT_READ_LIST=read_list.txt/" config.txt
        sed -i "s/GENOME_SIZE=/GENOME_SIZE="$genome_size"/" config.txt
        sed -i "s/THREADS=4/THREADS="$threads"/" config.txt
        ${netcat_dir}/necat.pl bridge config.txt
        cp necat/6-bridge_contigs/polished_contigs.fasta ${nextdenovo_dir}/${ID}/assembly_"$i".fasta
        rm -rf necat config.txt read_list.txt
    done

    for i in 05 11 17 23; do
        echo ${subsample_dir}/${ID}/sample_"$i".fastq > input.fofn
        cp "$NextDenovo_dir"/doc/run.cfg nextdenovo_run.cfg
        sed -i "s/genome_size = 1g/genome_size = "$genome_size"/" nextdenovo_run.cfg
        sed -i "s/parallel_jobs = 20/parallel_jobs = 1/" nextdenovo_run.cfg
        sed -i "s/read_type = clr/read_type = ont/" nextdenovo_run.cfg
        sed -i "s/pa_correction = 3/pa_correction = 1/" nextdenovo_run.cfg
        sed -i "s/correction_options = -p 15/correction_options = -p "$threads"/" nextdenovo_run.cfg
        sed -i "s/-t 8/-t "$threads"/" nextdenovo_run.cfg
        "$NextDenovo_dir"/nextDenovo nextdenovo_run.cfg
        cp 01_rundir/03.ctg_graph/nd.asm.fasta nextdenovo_temp.fasta
        rm -rf 01_rundir nextdenovo_run.cfg input.fofn
        echo ${subsample_dir}/${ID}/sample_"$i".fastq > lgs.fofn
        cat "$NextPolish_dir"/doc/run.cfg | grep -v "sgs" | grep -v "hifi" > nextpolish_run.cfg
        sed -i "s/parallel_jobs = 6/parallel_jobs = 1/" nextpolish_run.cfg
        sed -i "s/multithread_jobs = 5/multithread_jobs = "$threads"/" nextpolish_run.cfg
        sed -i "s|genome = ./raw.genome.fasta|genome = nextdenovo_temp.fasta|" nextpolish_run.cfg
        sed -i "s|-x map-ont|-x map-ont -t "$threads"|" nextpolish_run.cfg
        "$NextPolish_dir"/nextPolish nextpolish_run.cfg
        cp 01_rundir/genome.nextpolish.fasta ${nextdenovo_dir}/${ID}/assembly_"$i".fasta
        rm -rf 01_rundir pid*.log.info nextpolish_run.cfg lgs.fofn nextdenovo_temp.fasta
    done

    for i in 06 12 18 24; do
        raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly ${nextdenovo_dir}/${ID}/assembly_"$i".gfa ${subsample_dir}/${ID}/sample_"$i".fastq > ${nextdenovo_dir}/${ID}/assembly_"$i".fasta
    done
    } &
done

wait