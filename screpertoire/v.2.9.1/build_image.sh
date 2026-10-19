bsub -G compute-oncology -q oncology-interactive -Is -a 'docker_build(myhoang04/screpertoire:v2.9.1)' -- --tag myhoang04/screpertoire:v2.9.1 .
