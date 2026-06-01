cd genomics/

bsub -G compute-oncology -q general-interactive -Is -a 'docker_build(myhoang04/genomics)' -- --tag myhoang04/genomics:1.3 --tag myhoang04/genomics:latest .
