cd $(dirname $0)
bsub -G compute-oncology -q oncology-interactive -Is -a 'docker_build(myhoang04/splice2neo:v0.6.13_gencodemm10_ens95)' -- --tag myhoang04/splice2neo:v0.6.13_gencodemm10_ens95 .
