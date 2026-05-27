cd docker_image_builds/agfusion/agfusion_ens115
bsub -G compute-oncology -q general-interactive -Is -a 'docker_build(myhoang04/agfusion)' -- --tag myhoang04/agfusion:agfusion_ens115 .
