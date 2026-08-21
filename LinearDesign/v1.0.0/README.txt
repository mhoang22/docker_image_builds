cd docker_image_builds/lineardesign/v1.0.0
bsub -G compute-oncology -q general-interactive -Is -a 'docker_build(myhoang04/lineardesign)' -- --tag myhoang04/lineardesign:v1.0.0 .
