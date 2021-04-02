num_split=100 #split to sub files
num_sample=1000000 #random draw times
total_n=1732036 #snps
end=$((num_split-1))

for i in $(seq 0 $end)
#for i in {1..1}
do
    echo $i
    file='CHR_SGS_'$i'.sh'
    echo 'source ~/miniconda3/bin/activate' > $file
    echo 'Rscript SGS_HG_local.R '$i' '$num_split' '$num_sample' '$total_n >> $file
    chmod +x $file
    run_script $file

    #run_file='run_CHR_SGS_'$i'.sh'
    #echo 'EMAIL=`whoami`"@auburn.edu";' > $run_file
    #echo 'CWD=`pwd`;' >> $run_file
    #echo 'qsub -q general -N CHR_'$i' -j oe -e '$i'.error -l nodes=1:ppn=8,mem=32GB,walltime=36:00:00 -m be -M $EMAIL -d $CWD -V CHR_SGS_'$i'.sh' >> $run_file
    #echo 'exit 0' >> $run_file
done
