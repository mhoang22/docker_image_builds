bsub -G compute-oncology -q oncology-interactive -Is -a 'docker_build(myhoang04/diann:v2.5.1)' -- --tag myhoang04/diann:v2.5.1 .
