## Init objects to vdjtools scripts

#############
me=$(whoami)
application_files=/Users/$me/Dropbox/aa
files=/Users/$me/Dropbox/aa/

vdj=$application_files/applications/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjdb=$application_files/applications/vdjdb-1.1.5/vdjdb-1.1.5.jar
vdjdb_new=$application_files/data/selected_tcrb/databases/vdjdb_new


############# Init

cd $files/data/unselected_tcrb/Bethesda/; bethesda=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/Cleveland/; cleveland=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/ctrl_cd4/; healthy_cd4=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/ctrl_cd8/; healthy_cd8=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/Helsinki/; helsinki=$(ls -d "$PWD"/*); helsinki_pb=$(ls -d "$PWD"/*pb*);
cd $files/data/unselected_tcrb/Helsinki_itp/; helsinki_itp=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/Helsinki_bmf/; helsinki_bmf=$(ls -d "$PWD"/*); helsinki_bmf_pb=$(ls -d "$PWD"/*pb*);
cd $files/data/unselected_tcrb/Helsinki_RA/; helsinki_ra=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/Melbourne/; melbourne=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/Melbourne_bmfu/; melbourne_bmfu=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/Melbourne_ibmf/; melbourne_ibmf=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/Melbourne_mds/; melbourne_mds=$(ls -d "$PWD"/*);
cd $files/data/unselected_tcrb/Melbourne_pnh/; melbourne_pnh=$(ls -d "$PWD"/*);

###############
cd $files
clear
###############

## Preprocess; convert into more human readable format
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $bethesda data/unselected_tcrb/Bethesda/vdj/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $healthy_cd4 data/unselected_tcrb/Healthy/ctrl_cd4/vdj/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $healthy_cd8 data/unselected_tcrb/Healthy/ctrl_ctl/vdj/
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $helsinki_bmf data/unselected_tcrb/Helsinki_bmf/vdj/

## Remove unproductive clonotypes
java -Xmx4G -jar $vdj FilterNonFunctional $bethesda data/unselected_tcrb/Bethesda/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $cleveland data/unselected_tcrb/Cleveland/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $healthy_cd4 data/unselected_tcrb/Healthy/ctrl_cd4/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $healthy_cd8 data/unselected_tcrb/Healthy/ctrl_ctl/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $helsinki data/unselected_tcrb/Helsinki/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $helsinki_bmf data/unselected_tcrb/Helsinki_bmf/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $helsinki_ra data/unselected_tcrb/Helsinki_RA/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $melbourne data/unselected_tcrb/Melbourne/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $melbourne_bmfu data/unselected_tcrb/Melbourne_bmfu/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $melbourne_ibmf data/unselected_tcrb/Melbourne_ibmf/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $melbourne_mds data/unselected_tcrb/Melbourne_mds/filtered/;
java -Xmx4G -jar $vdj FilterNonFunctional $melbourne_pnh data/unselected_tcrb/Melbourne_pnh/filtered/;

## Downsample CD8-sorted samples to 10k
java -Xmx4G -jar $vdj DownSample --size 10000 $bethesda            data/unselected_tcrb/resampled/Bethesda/
java -Xmx4G -jar $vdj DownSample --size 10000 $healthy_cd4         data/unselected_tcrb/resampled/Healthy/ctrl_cd4/
java -Xmx4G -jar $vdj DownSample --size 10000 $healthy_cd8         data/unselected_tcrb/resampled/Healthy/ctrl_cd8/
java -Xmx4G -jar $vdj DownSample --size 10000 $helsinki_pb         data/unselected_tcrb/resampled/Helsinki/
java -Xmx4G -jar $vdj DownSample --size 10000 $helsinki_bmf_pb_dg  data/unselected_tcrb/resampled/Helsinki_bmf/
java -Xmx4G -jar $vdj DownSample --size 10000 $helsinki_ra         data/unselected_tcrb/resampled/Helsinki_RA/
java -Xmx4G -jar $vdj DownSample --size 10000 $helsinki_itp        data/unselected_tcrb/resampled/Helsinki_itp/

## Downsample MNC-sorted samples to 20k
java -Xmx4G -jar $vdj DownSample --size 20000 $cleveland      data/unselected_tcrb/resampled/Cleveland/
java -Xmx4G -jar $vdj DownSample --size 20000 $melbourne      data/unselected_tcrb/resampled/Melbourne/
java -Xmx4G -jar $vdj DownSample --size 20000 $melbourne_bmfu data/unselected_tcrb/resampled/Melbourne_bmfu/
java -Xmx4G -jar $vdj DownSample --size 20000 $melbourne_ibmf data/unselected_tcrb/resampled/Melbourne_ibmf/
java -Xmx4G -jar $vdj DownSample --size 20000 $melbourne_mds  data/unselected_tcrb/resampled/Melbourne_mds/
java -Xmx4G -jar $vdj DownSample --size 20000 $melbourne_pnh  data/unselected_tcrb/resampled/Melbourne_pnh/
java -Xmx4G -jar $vdj DownSample --size 20000 $emerson        data/unselected_tcrb/resampled/emerson/

## CalcDiversityStats
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa - $bethesda $cleveland $healthy_cd4 $healthy_cd8 $helsinki $itp $helsinki_bmf $helsinki_ra $melbourne $melbourne_bmfu $melbourne_ibmf $melbourne_mds $melbourne_pnh $files/tcrb_data/results/diversity/
